"""Plotting utilities for ct_validation enrichment results."""

from __future__ import annotations

from typing import TYPE_CHECKING, Literal

if TYPE_CHECKING:
    import matplotlib.pyplot as plt
    import pandas as pd

_METRIC_COLS = {
    "rr": ("rr", "rr_ci_lower", "rr_ci_upper", "Risk Ratio"),
    "or": ("or", "or_ci_lower", "or_ci_upper", "Odds Ratio"),
}


def forest_plot(
    data: pd.DataFrame,
    label_col: str = "phase_label",
    phase: str | None = None,
    metric: Literal["rr", "or"] = "rr",
    title: str = "",
    figsize: tuple[float, float] = (5, 3),
    ax: plt.Axes | None = None,
    *,
    show_counts: bool = True,
    sort_order: list[str] | None = None,
    color_col: str | None = None,
) -> tuple[plt.Figure, plt.Axes]:
    """Horizontal forest plot of effect sizes with 95% confidence intervals.

    Parameters
    ----------
    data : pd.DataFrame
        Enrichment results from ``validate()``. Must contain the columns
        for the chosen *metric* (e.g. ``rr``, ``rr_ci_lower``,
        ``rr_ci_upper`` for risk ratio) and the column named by
        *label_col*. If *show_counts* is True, also needs ``x_yes``
        and ``n_yes``.
    label_col : str
        Column used for y-axis labels (default ``"phase_label"``).
    phase : str or None
        If given, filter ``data`` to rows where ``phase_label == phase``.
    metric : ``"rr"`` or ``"or"``
        Which effect size to plot (default ``"rr"``).
    title : str
        Axes title.
    figsize : tuple
        Figure size when creating a new figure (ignored if *ax* is given).
    ax : matplotlib Axes or None
        Axes to draw on. If None a new figure is created.
    show_counts : bool
        Annotate each point with ``x_yes/n_yes`` counts (default True).
    sort_order : list of str or None
        Explicit label order for the y-axis. If None (default), rows are
        sorted by the metric ascending.
    color_col : str or None
        Column whose values split the rows into colour groups, drawn with
        a legend. Row order is unaffected.

    Returns
    -------
    (fig, ax)
    """
    try:
        import matplotlib.pyplot as plt  # noqa: PLC0415
    except ModuleNotFoundError:
        msg = (
            "matplotlib is required for plotting. Install it with: pip install ct-validation[plot]"
        )
        raise ModuleNotFoundError(msg) from None

    col, ci_lo, ci_hi, xlabel = _METRIC_COLS[metric]

    df = data.copy()
    if phase:
        df = df[df["phase_label"] == phase]
    if sort_order is not None:
        df = df.set_index(label_col).loc[sort_order].reset_index()
    else:
        df = df.sort_values(col, ascending=True).reset_index(drop=True)

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.figure

    y = range(len(df))
    # y positions come from the sorted frame's index, so colour groups may interleave.
    # dropna=False: a dropped row loses its marker while its tick label and count still
    # draw, which reads as an off-scale point.
    groups = [(None, df)] if color_col is None else df.groupby(color_col, sort=False, dropna=False)
    for name, part in groups:
        ax.errorbar(
            part[col],
            part.index,
            xerr=[part[col] - part[ci_lo], part[ci_hi] - part[col]],
            fmt="o",
            label=None if name is None else str(name),
        )
    if color_col is not None:
        ax.legend(frameon=False, fontsize=9)
    ax.set_yticks(list(y))
    ax.set_yticklabels(df[label_col].tolist())
    ax.axvline(1.0, color="lightgray", linestyle="--", linewidth=0.8)
    ax.set_xlabel(xlabel)
    if title:
        ax.set_title(title)

    if show_counts:
        for i, (_, row) in enumerate(df.iterrows()):
            ax.annotate(
                f"{row['x_yes']:,}/{row['n_yes']:,}",
                (1.02, i),
                xycoords=("axes fraction", "data"),
                fontsize=10,
                va="center",
                ha="left",
                annotation_clip=False,
            )

    return fig, ax
