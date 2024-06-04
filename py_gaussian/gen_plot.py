import matplotlib as plt
import pandas as pd
import numpy as np
from scipy.optimize import minimize
from scipy import interpolate
from exceptions import CouldNotPlotSmoothProfile
import textalloc as ta


def plot_smooth_profile(zi_s, energies, ax, fig, label, color=None):
    """
    Plot a smooth reaction profile by spline interpolation and finding the
    stationary points. This will not afford the correct number of stationary
    points for some energy arrays, so raise an exception if it fails

    ---------------------------------------------------------------------------
    Arguments:
        zi_s (np.ndarray): Estimate of reaction coordinate points

        energies (list(autode.plotting.Energy)): len(energies) = len(zi_s)

        ax (matplotlib.axes.Axes):
    """

    # Minimise a set of spline points so the stationary points have y values
    # given in the energies array
    energies_arr = np.array([energy for energy in energies], dtype="f")
    result = minimize(
        error_on_stationary_points,
        x0=energies_arr,
        args=(energies_arr,),
        method="BFGS",
        tol=0.1,
    )

    # Use the optimised values to construct a spline function that will be
    # plotted
    optimised_spline = interpolate.CubicSpline(
        zi_s, result.x, bc_type="clamped"
    )

    # Create more zi values from slightly before the minimum to slightly after
    # the maximum
    fine_zi_s = np.linspace(min(zi_s) - 0.2, max(zi_s) + 0.2, num=500)

    # The new zi values are the stationary points of the optimised function
    zi_s = get_stationary_points(fine_zi_s, optimised_spline.derivative())

    if len(zi_s) != len(energies):
        raise CouldNotPlotSmoothProfile

    # Plot the function
    if color:
        ax.plot(fine_zi_s, optimised_spline(fine_zi_s), label=label, color=color)
    else:
        ax.plot(fine_zi_s, optimised_spline(fine_zi_s), label=label)
    #ax.scatter(zi_s, optimised_spline(zi_s), zorder=10)

    # Annotate the plot with the relative energies
    max_delta = max(energies) - min(energies)
    
    txt_lst = []
    for i, energy in enumerate(optimised_spline(zi_s)):
        # Shift the minima labels (even points) below the point and the
        # transition state labels above the point
        #shift = -0.07 * max_delta if i % 2 == 0 else 0.03 * max_delta

        """ax.annotate(
            f"{energy:.1f}",
            (zi_s[i], energy + shift),
            fontsize=12,
            ha="center",
        )"""
        txt_lst.append(f"{energy:.1f}")
    """ta.allocate_text(fig,ax,zi_s,optimised_spline(zi_s),
                txt_lst,
                x_scatter=zi_s, y_scatter=optimised_spline(zi_s),
                x_lines=[fine_zi_s],y_lines=[optimised_spline(fine_zi_s)],
                draw_lines=False,
                textsize=10)"""


    return None


def error_on_stationary_points(x, energies):
    """
    Calculate the difference between the stationary points of an interpolated
    function and those observed (given in the energies array). Example::

          |      .
        E |.   / |        The points indicate the true stationary points
          | |_/  |.
          |_____________
                zi

    ---------------------------------------------------------------------------
    Arguments:
        x (np.ndarray): Points that will be splined that generate stationary
                        points that ≈ energies

        energies (np.ndarray): Observed stationary points

    Returns:
        (float): A measure of the error
    """
    # Generate a list of reaction coordinate points - arbitrary units so
    # integers are fine
    zi_s = np.array(range(len(x)))

    # Spline the energies to get a function that has stationary points
    spline = interpolate.CubicSpline(zi_s, x, bc_type="clamped")

    # Calculate the energy values at the stationary points of the function with
    # a fine-ish spacing that extrapolates
    # slightly
    fine_zi_s = np.linspace(min(zi_s) - 0.2, max(zi_s) + 0.2, num=500)
    stationary_points = get_stationary_points(
        xs=fine_zi_s, dydx=spline.derivative()
    )

    if len(stationary_points) != len(energies):
        # TODO make this smooth somehow
        # Energy penalty for not having the required number of
        return 10 * np.abs(len(energies) - len(stationary_points))

    energies_at_stationary_points = [spline(zi) for zi in stationary_points]

    # Return the error as the sum squared difference between the required and
    # the observed stationary point energies
    energy_difference = energies - np.array(energies_at_stationary_points)

    return np.sum(np.square(energy_difference))


def get_stationary_points(xs, dydx):
    """
    Compute the productive of the derivative at points x(i-1) and x(i) which
    is negative if there is a point x(k)
    between x(i-1) and x(i) that has dy/dx|x(k) = 0

    ---------------------------------------------------------------------------
    Arguments:
         xs (np.ndarray):

         dydx (function):
    """
    stationary_points = []

    for i in range(1, len(xs) - 1):
        if dydx(xs[i - 1]) * dydx(xs[i]) < 0:
            stationary_points.append(xs[i])

    return stationary_points