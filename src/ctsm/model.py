"""Multilayer perceptron classifier and its training loop (decision D-08).

Architecture: input dropout -> [Linear -> BatchNorm1d -> ReLU -> Dropout] x n_layers -> Linear.
Hidden widths taper by half per layer (width, width/2, width/4; minimum 16). Every layer type is
supported exactly by shap.DeepExplainer, and each module instance is used once.
Outputs are logits: 1 (binary, log-odds of tumour) or C (multiclass).
Loss: class-balanced cross-entropy (each class contributes equally regardless of size).
Early stopping on the validation balanced log loss; the best epoch's weights are restored.
"""
from __future__ import annotations

import copy

import numpy as np
import torch
from torch import nn


class MLP(nn.Module):
    def __init__(self, n_in: int, n_out: int, n_layers: int, width: int, dropout: float, input_dropout: float):
        super().__init__()
        layers: list[nn.Module] = [nn.Dropout(input_dropout)]
        d = n_in
        for i in range(n_layers):
            w = max(16, int(width) // (2 ** i))
            layers += [nn.Linear(d, w), nn.BatchNorm1d(w), nn.ReLU(), nn.Dropout(dropout)]
            d = w
        layers.append(nn.Linear(d, n_out))
        self.net = nn.Sequential(*layers)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return self.net(x)


def class_weights(y: np.ndarray, n_classes: int) -> np.ndarray:
    """n / (C_present * n_c): every class present in y carries equal total weight."""
    counts = np.bincount(y, minlength=n_classes).astype(np.float64)
    present = counts > 0
    w = np.zeros(n_classes)
    w[present] = len(y) / (present.sum() * counts[present])
    return w


def _loss(logits: torch.Tensor, y: torch.Tensor, sw: torch.Tensor, binary: bool) -> torch.Tensor:
    if binary:
        per = nn.functional.binary_cross_entropy_with_logits(logits[:, 0], y.float(), reduction="none")
    else:
        per = nn.functional.cross_entropy(logits, y, reduction="none")
    return (per * sw).sum() / sw.sum()


def balanced_log_loss(logits: np.ndarray, y: np.ndarray, binary: bool) -> float:
    """Mean over classes present in y of the class-mean cross-entropy (a proper scoring rule
    that is insensitive to class imbalance)."""
    t = torch.as_tensor(logits, dtype=torch.float64)
    yt = torch.as_tensor(y)
    if binary:
        per = nn.functional.binary_cross_entropy_with_logits(t[:, 0], yt.double(), reduction="none").numpy()
    else:
        per = nn.functional.cross_entropy(t, yt, reduction="none").numpy()
    return float(np.mean([per[y == c].mean() for c in np.unique(y)]))


def predict_logits(model: nn.Module, X: np.ndarray, batch: int = 4096) -> np.ndarray:
    model.eval()
    out = []
    with torch.no_grad():
        for a in range(0, len(X), batch):
            out.append(model(torch.from_numpy(np.ascontiguousarray(X[a:a + batch]))).numpy())
    return np.concatenate(out, axis=0)


def logits_to_proba(logits: np.ndarray, binary: bool) -> np.ndarray:
    """Binary: P(tumour) as an (n, 2) array [P(normal), P(tumour)]; multiclass: softmax."""
    if binary:
        p1 = 1.0 / (1.0 + np.exp(-logits[:, 0].astype(np.float64)))
        return np.column_stack([1 - p1, p1])
    z = logits.astype(np.float64)
    z = z - z.max(axis=1, keepdims=True)
    e = np.exp(z)
    return e / e.sum(axis=1, keepdims=True)


def train_model(X_tr: np.ndarray, y_tr: np.ndarray, params: dict, n_classes: int, seed: int, *,
                X_val: np.ndarray | None = None, y_val: np.ndarray | None = None,
                max_epochs: int = 200, patience: int = 20, fixed_epochs: int | None = None):
    """Train one MLP. With a validation set: early stopping, returns the best-epoch weights.
    With fixed_epochs: trains exactly that many epochs (final refit on the full outer-train).
    Returns (model, best_epoch [1-based], best_val_loss or nan, history)."""
    binary = n_classes == 2
    torch.manual_seed(seed)
    rng = np.random.default_rng(seed)
    model = MLP(X_tr.shape[1], 1 if binary else n_classes, params["n_layers"], params["width"],
                params["dropout"], params["input_dropout"])
    opt = torch.optim.AdamW(model.parameters(), lr=params["lr"], weight_decay=params["weight_decay"])
    Xt = torch.from_numpy(np.ascontiguousarray(X_tr, dtype=np.float32))
    yt = torch.from_numpy(y_tr.astype(np.int64))
    sw = torch.from_numpy(class_weights(y_tr, n_classes)[y_tr].astype(np.float32))
    bs = int(params["batch_size"])
    n = len(y_tr)
    epochs = fixed_epochs if fixed_epochs is not None else max_epochs
    best = (np.inf, 0, None)
    history = []
    for ep in range(1, epochs + 1):
        model.train()
        perm = torch.from_numpy(rng.permutation(n))
        for a in range(0, n, bs):
            idx = perm[a:a + bs]
            if len(idx) < 2:                 # BatchNorm needs >= 2 samples
                continue
            opt.zero_grad(set_to_none=True)
            loss = _loss(model(Xt[idx]), yt[idx], sw[idx], binary)
            loss.backward()
            opt.step()
        if X_val is None:
            continue
        v = balanced_log_loss(predict_logits(model, X_val), y_val, binary)
        v = v if np.isfinite(v) else np.inf              # diverged: never "best"
        history.append(v)
        if v < best[0] - 1e-6:
            best = (v, ep, copy.deepcopy(model.state_dict()))
        elif ep - best[1] >= patience:
            break
    if X_val is not None:
        if best[2] is None:                              # no finite validation loss at all
            return model, len(history), np.inf, history
        model.load_state_dict(best[2])
        return model, best[1], best[0], history
    return model, epochs, float("nan"), history
