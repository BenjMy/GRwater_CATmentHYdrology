"""
agramon_diagnostics.py - forcing diagnostics figures for the Agramon pipeline
============================================================================
Moved unchanged from Agramon_withLAI_withETp.py (REFACTORING.md 2.6); the main script
re-exports the name, so existing imports keep working.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def save_forcing_diagnostics(
    valid_times,
    net_flux_1d,
    rain_nodes=None,
    ETp_nodes=None,
    net_flux_2d=None,
    outdir="forcing_plots",
):
    """
    Save forcing diagnostic figures.

    Parameters
    ----------
    valid_times : array-like
        Time axis.
    net_flux_1d : ndarray (nt,)
        Mean net flux [m/s].
    rain_nodes : ndarray (nt, nnodes), optional
        Rainfall forcing [m/s].
    ETp_nodes : ndarray (nt, nnodes), optional
        ETp forcing [m/s].
    net_flux_2d : ndarray (nt, nnodes), optional
        Spatially distributed net flux [m/s].
    outdir : str or Path
        Output directory.
    """

    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # --------------------------------------------------
    # Net flux timeseries
    # --------------------------------------------------
    net_mmday = net_flux_1d * 86400 * 1000

    fig, ax = plt.subplots(figsize=(15, 4))

    ax.plot(
        valid_times,
        net_mmday,
        color="k",
        lw=1.5,
        label="Rain - ETp",
    )

    ax.axhline(0, color="0.5", ls="--", lw=0.8)

    ax.fill_between(
        valid_times,
        net_mmday,
        0,
        where=(net_mmday > 0),
        color="royalblue",
        alpha=0.4,
        label="Net input",
    )

    ax.fill_between(
        valid_times,
        net_mmday,
        0,
        where=(net_mmday < 0),
        color="firebrick",
        alpha=0.4,
        label="Net ET loss",
    )

    ax.set_title("Net Atmospheric Flux")
    ax.set_ylabel("mm/day")
    ax.grid(alpha=0.3)
    ax.legend()

    fig.tight_layout()
    fig.savefig(outdir / "net_flux_timeseries.png",
                dpi=300, bbox_inches="tight")
    plt.close(fig)

    # --------------------------------------------------
    # Components
    # --------------------------------------------------
    if rain_nodes is not None and ETp_nodes is not None:

        mean_rain = rain_nodes.mean(axis=1) * 86400 * 1000
        mean_etp = ETp_nodes.mean(axis=1) * 86400 * 1000

        fig, ax = plt.subplots(figsize=(15, 5))

        ax.plot(valid_times, mean_rain,
                color="blue", label="Rain")

        ax.plot(valid_times, mean_etp,
                color="orange", label="ETp")

        ax.plot(valid_times, net_mmday,
                color="black", lw=2,
                label="Rain - ETp")

        ax.axhline(0, color="0.5", ls="--")

        ax.set_title("Atmospheric Forcing Components")
        ax.set_ylabel("mm/day")
        ax.grid(alpha=0.3)
        ax.legend()

        fig.tight_layout()
        fig.savefig(outdir / "forcing_components.png",
                    dpi=300, bbox_inches="tight")
        plt.close(fig)

    # --------------------------------------------------
    # Spatial variability
    # --------------------------------------------------
    if net_flux_2d is not None:

        p10 = np.percentile(net_flux_2d, 10, axis=1) * 86400 * 1000
        p50 = np.percentile(net_flux_2d, 50, axis=1) * 86400 * 1000
        p90 = np.percentile(net_flux_2d, 90, axis=1) * 86400 * 1000

        fig, ax = plt.subplots(figsize=(15, 4))

        ax.plot(valid_times, p50,
                color="k", label="Median")

        ax.fill_between(
            valid_times,
            p10,
            p90,
            alpha=0.3,
            label="10–90 %",
        )

        ax.set_title("Spatial Distribution of Net Flux")
        ax.set_ylabel("mm/day")
        ax.grid(alpha=0.3)
        ax.legend()

        fig.tight_layout()
        fig.savefig(
            outdir / "net_flux_spatial_variability.png",
            dpi=300,
            bbox_inches="tight",
        )
        plt.close(fig)

        # --------------------------------------------------
        # Annual water balance
        # --------------------------------------------------
        if rain_nodes is not None and ETp_nodes is not None:

            import pandas as pd

            # Mean forcing over surface nodes
            rain_mmday = rain_nodes.mean(axis=1) * 86400 * 1000
            etp_mmday  = ETp_nodes.mean(axis=1) * 86400 * 1000
            net_mmday  = rain_mmday - etp_mmday

            df = pd.DataFrame(
                {
                    "rain": rain_mmday,
                    "etp": etp_mmday,
                    "net": net_mmday,
                },
                index=pd.to_datetime(valid_times),
            )

            # Annual totals [mm/year]
            annual = df.resample("Y").sum()
            annual.index = annual.index.year

            fig, ax = plt.subplots(figsize=(10, 5))

            x = np.arange(len(annual))
            w = 0.28

            ax.bar(
                x - w,
                annual["rain"],
                width=w,
                color="royalblue",
                label="Rain",
            )

            ax.bar(
                x,
                annual["etp"],
                width=w,
                color="darkorange",
                label="ETp",
            )

            ax.bar(
                x + w,
                annual["net"],
                width=w,
                color="forestgreen",
                label="Net",
            )

            ax.axhline(0, color="k", lw=0.8)

            ax.set_xticks(x)
            ax.set_xticklabels(annual.index.astype(str))

            ax.set_ylabel("Annual total [mm]")
            ax.set_title("Annual Water Balance")
            ax.legend()

            fig.tight_layout()

            fig.savefig(
                outdir / "annual_water_balance.png",
                dpi=300,
                bbox_inches="tight",
            )

            plt.close(fig)

            years = np.unique(pd.to_datetime(valid_times).year)

            ncols = 3
            nrows = int(np.ceil(len(years) / ncols))

            fig, axes = plt.subplots(
                nrows,
                ncols,
                figsize=(4*ncols, 3*nrows),
                sharex=True,
                sharey=True,
            )

            axes = np.atleast_1d(axes).ravel()

            for ax, year in zip(axes, years):

                mask = pd.to_datetime(valid_times).year == year

                ax.hist(
                    net_mmday[mask],
                    bins=40,
                    color="steelblue",
                    alpha=0.8,
                )

                ax.axvline(0, color="k", ls="--", lw=1)
                ax.set_title(str(year))

            for ax in axes[len(years):]:
                ax.set_visible(False)

            fig.suptitle("Distribution of Daily Net Flux by Year")
            fig.tight_layout()

            fig.savefig(
                outdir / "net_flux_histogram_by_year.png",
                dpi=300,
                bbox_inches="tight",
            )

            plt.close(fig)

            # ==================================================
            # Monthly climatology
            # ==================================================

            monthly = (
                df.groupby(df.index.month)
                  .sum()
            )

            month_labels = [
                "Jan", "Feb", "Mar", "Apr", "May", "Jun",
                "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"
            ]

            fig, ax = plt.subplots(figsize=(12, 5))

            x = np.arange(12)
            w = 0.28

            ax.bar(
                x - w,
                monthly["rain"],
                width=w,
                color="royalblue",
                label="Rain",
            )

            ax.bar(
                x,
                monthly["etp"],
                width=w,
                color="darkorange",
                label="ETp",
            )

            ax.bar(
                x + w,
                monthly["net"],
                width=w,
                color="forestgreen",
                label="Net",
            )

            ax.axhline(0, color="k", lw=0.8)

            ax.set_xticks(x)
            ax.set_xticklabels(month_labels)

            ax.set_ylabel("Accumulated forcing [mm]")
            ax.set_title("Monthly Water Balance (all years combined)")
            ax.legend()

            fig.tight_layout()

            fig.savefig(
                outdir / "monthly_water_balance.png",
                dpi=300,
                bbox_inches="tight",
            )

            plt.close(fig)

            fig, ax = plt.subplots(figsize=(12, 5))

            data = [
                df.loc[df.index.month == m, "net"].values
                for m in range(1, 13)
            ]

            ax.boxplot(
                data,
                labels=month_labels,
                showfliers=False,
            )

            ax.axhline(0, color="k", ls="--")

            ax.set_ylabel("Net flux [mm/day]")
            ax.set_title("Monthly Distribution of Net Flux")

            fig.tight_layout()

            fig.savefig(
                outdir / "monthly_net_flux_boxplot.png",
                dpi=300,
                bbox_inches="tight",
            )

            plt.close(fig)
