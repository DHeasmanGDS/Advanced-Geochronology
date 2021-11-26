"""
Author                      : Drew Heasman
Date(last updated)          : 05 June 2021
Description : Advanced Geochronology code. This code is designed to be used with the Jupyter Notebooks created for Geochron calcualtions given by Bruce Eglington and Camille Partin for the Advanced Geochronology class 
            at the University of Saskatchewan.
"""

import pandas as pd
import numpy as np
from bokeh.plotting import figure
from bokeh.models import ColumnDataSource, LabelSet
from scipy.interpolate import CubicSpline

def get_figure(title, x_label, y_label, x_range, y_range):
    """Generates a blank figure for Concordia Diagrams.

    Parameters
    ----------
    title : string
        The title of the figure.
    x_label : string
        Label of the x-axis.
    y_label : string
        Label of the y-axis.
    x_range : list
        Range for the x-axis, a minimum and maximum.
    y_range : list
        Range for the y-axis, a minimum and maximum.

    Returns
    -------
    figure
        A template figure.

    Raises
    ------
    None
    """
    fig = figure(title=title,
           x_axis_label=x_label,
           y_axis_label=y_label,
           width=700,
           plot_height=500,
           toolbar_location="left",
           toolbar_sticky=False,
           y_range=y_range,
           x_range=x_range,
           tools=['pan,box_zoom,crosshair,hover,reset,save']
           )
    fig.toolbar.logo = None
    fig.title.align = "center"
    fig.title.text_font_size = "25px"
    return fig

def plot_wetherill(fig, df, color, legend_label, width):
    """Plots the Wetherill Concordia line on a figure.

    Parameters
    ----------
    fig : figure
        The figure you want to plot the diagram on.
    df : DataFrame
        A dataframe of the data, must have columns "ratio_207Pb_235U" and "ratio_206Pb_238U"
    color : string
        The color of the plotting line.
    legend_label : string
        What to label that line.
    width : int
        The thickness of the line.

    Returns
    -------
    figure 
        A figure with a plot of the Wetherill Concordia Line.

    Raises
    ------
    None.
    """
    fig.line(df["ratio_207Pb_235U"], df["ratio_206Pb_238U"],\
            color=color, legend_label=legend_label, line_width=width)
    return fig

def plot_ages_wetherill(fig, ages):
    """Plots points of ages on a graph along the Wetherill Concordia Line.
    
    Parameters
    ----------
    fig : figure
        The figure you want to plot the diagram on.
    ages : list
        A list of ages you want to add points to on a Wetherill plot.
        
    Returns
    -------
    figure 
        A figure with age points on the Wetherill line.

    Raises
    ------
    None.
    """
    source = ColumnDataSource(ages)
    labels = LabelSet(x='ratio_207Pb_235U', y='ratio_206Pb_238U', text='age', level='overlay',\
              x_offset=10, y_offset=-10, source=source, render_mode='canvas')
    fig.add_layout(labels)
    return fig

def plot_tera_wasserburg(fig, df, color, legend_label, width):
    """Plots the Tera-Wasserburg Concordia line on a figure.

    Parameters
    ----------
    fig : figure
        The figure you want to plot the diagram on.
    df : DataFrame
        A dataframe of the data, must have columns "ratio_207Pb_235U" and "ratio_206Pb_238U"
    color : string
        The color of the plotting line.
    legend_label : string
        What to label that line.
    width : int
        The thickness of the line.

    Returns
    -------
    figure 
        A figure with a Tera-Wasserburg Concordia plot.

    Raises
    ------
    None.
    """
    fig.line(df["ratio_238U_206Pb"], df["ratio_207Pb_206Pb"],\
            color=color, legend_label=legend_label, line_width=width)
    return fig

def plot_ages_tera_wasserburg(fig, ages):
    """Plots points of ages on a graph along the Tera-Wasserburg Line.
    
    Parameters
    ----------
    fig : figure
        The figure you want to plot the diagram on.
    ages : list
        A list of ages you want to add points to on a Tera-Wasserburg plot.
        
    Returns
    -------
    figure 
        A figure with age points on the Tera-Wasserburg line.

    Raises
    ------
    None.
    """
    source = ColumnDataSource(ages)
    labels = LabelSet(x='ratio_238U_206Pb', y='ratio_207Pb_206Pb', text='age', level='overlay',\
              x_offset=10, y_offset=-10, source=source, render_mode='canvas')
    fig.add_layout(labels)
    return fig


def get_U_ratio_data(lambda238=1.55125*10**-10, lambda235=9.8485*10**-10, age=np.append([0.001], [np.arange(1,4568)])):
    """Generates a DataFrame with various ratios for concordia plots.

    Parameters
    ----------
    lambda238 : float (Default value = 1.55125*10**-10)
        The decay constant for Uranium 238.
    lambda235 : float (Default value = 9.8485*10**-10)
        The decay constant for Uranium 235
    age : (Default value = np.append([0.001]), [np.arange(1,4568)])
        A list with a range, per year starting at 0.001, and then every integer between 1 and 4568.
        
    Returns
    -------
    pandas DataFrame
        A dataframe with uranium ration data for calculation of various concordia plots. 
    """
    age = age
    df = pd.DataFrame()
    df["age"] = age
    df["ratio_238U_206Pb"] = 1/(np.exp(lambda238*df.age*1000000)-1)
    df["ratio_207Pb_235U"] = np.exp(lambda235*df.age*1000000)-1
    df["ratio_206Pb_238U"] = np.exp(lambda238*df.age*1000000)-1
    df["ratio_207Pb_206Pb"] = (1/137.88)*(df.ratio_207Pb_235U/df.ratio_206Pb_238U)    
    return df

def calc_t2_daughter(initial,
                  t2_parent,
                  decay_constant,
                  t1_age,
                  t2_age=0):
    """Gives the present day (or t2) daughter composition, using the isochron equation.

    Parameters
    ----------
    initial : float
        The initial daughter ratio (or t1). 
    t2_parent : float
        The t2 parent ratio.
    decay_constant : float
        The decay constant of the parent.
    t1_age : float
        The initial (or t1) age.
    t2_age : float (Default value = 0)
        The t2 age, which defaults to present day

    Returns
    -------
    float
        The t2 daughter composition.
        
    Raises
    ------
    None.
    """
    return initial + (t2_parent * (np.exp(decay_constant*t1_age*1000000)\
                                  -np.exp(decay_constant*t2_age*1000000)))

def calc_initial(t2_daughter,
                  t2_parent,
                  decay_constant,
                  t1_age,
                  t2_age=0):
    """Gives the initial (or t1) daughter composition, using the isochron equation.
    Parameters
    ----------
    t2_daughter : float
        The t2 daughter ratio.
    t2_parent : float
        The t2 parent ratio.
    decay_constant : float
        The decay constant of the parent
    t1_age : float
        The initial (or t1) age.
    t2_age : float (Default value = 0)
        The t2 age, which defaults to present day.
    
    Returns
    -------
    float
        The initial (or t1) daughter composition
    
    Raises
    ------
    None.
    """
    return t2_daughter - (t2_parent * (np.exp(decay_constant*t1_age*1000000)\
                                  -np.exp(decay_constant*t2_age*1000000)))

def calc_age(est_t2_daughter,
             t2_parent,
             t2_daughter,
             decay_constant):
    """Returns an age, using the isochron equation.

    Parameters
    ----------
    est_initial : float
        The estimated initial daughter ratio.
    t2_parent : float
        The measured t2 parent ratio.
    present_day : float
        The measured t2 daughter ratio.
    decay_constant : float
        The decay constant of the parent.

    Returns
    -------
    float
        An apparent age in Ma.
        
    Raises
    ------
    None.
    """
    return np.round((np.log(((t2_daughter - est_t2_daughter)/t2_parent + 1)) / decay_constant) / 1000000, 1)

def convert_halflife(halflife):
    """Returns a decay constant given a half life.

    Parameters
    ----------
    halflife : float
        A half life of a given isotope.

    Returns
    -------
    float
        The decay constant of that isotope.
    
    Raises
    ------
    None.
    """
    return 0.693 / halflife / 1000000

def get_age_df(df, age_list):
    """Gets a DataFrame and filters it for the purpose of putting specific ages on the Concordia plots. Note the DataFrame must have a column named 'age'.

    Parameters
    ----------
    df : DataFrame
        A dataframe with a column 'age'.
    age_list : list
        A list of ages you want to place on a Concordia plot.

    Returns
    -------
    DataFrame 
        A filtered pandas DataFrame.

    Raises
    ------
    None.
    """
    boolean_series = df.age.isin(age_list)
    return df[boolean_series]

def get_PbPb_ratio_data(mu,
                        kappa,
                        t1=3700,
                        t2=0,
                        Pb206i=11.152,
                        Pb207i=12.988,
                        Pb208i=31.23):
    """Generates a DataFrame with various ratios for PbPb plots.

    Parameters
    ----------
    t1 : int (Default value = 3700)
        The max age for calculating the PbPb data.  The default is based on Stacey and Kramers (1975)
    t2 : int (Default value = 0)
        The min age for calculating the PbPb data.
    mu : float
        U238/204Pb ratio.  
    Pb206i : float (Default value = 11.152)
        Pb206/Pb204 initial ratio. Default is based on Stacey and Kramers (1975)
    Pb207i : float (Default value = 12.988)
        Pb207/Pb204 initial ratio. Default is based on Stacey and Kramers (1975)
    Pb208i : float (Default value = 31.23)
        Pb208/Pb204 initial ratio. Default is based on Stacey and Krameres (1975)
    kappa : float
        232Th/238U ratio, which depends on the source.  A good number for depleted mantle is 2.5, ocean island basalts is 3.75, lower mantle is 4, and continental crust is 5. (Paul, White, and Turcotte 2003)

    Returns
    -------
    DataFrame
        A dataframe useful for plotting Pb/Pb curves.

    Raises
    ------
    None.
    """

    decay_const_238 = 1.55125*10**-10   # 1975
    decay_const_235 = 9.8485*10**-10    # 1975
    decay_const_232 = 4.9745*10**-11    # 1975
    U238U235 = 137.88
    d=[]
    while t2 < t1:
        
        d.append({'t1': t1,
            't2': t2,
            '238U/204Pb': mu, 
            '206Pb/204Pb': calc_t2_daughter(Pb206i, mu, decay_const_238, t1, t2),
            '207Pb/204Pb': calc_t2_daughter(Pb207i, mu/U238U235, decay_const_235, t1, t2),
            '208Pb/204Pb': calc_t2_daughter(Pb208i, kappa * mu, decay_const_232, t1, t2)
                 })
        t2 = t2 + 1
    return pd.DataFrame(d).sort_values(by='206Pb/204Pb')

def plot_uranogenicPb_curve(fig,df,color,label):
    """Plot a uranogenic Pb/Pb curve.  206Pb/204Pb vs. 207Pb/204Pb.
    Parameters
    ----------
    fig : figure
        The figure you want to plot the diagram on.
    df : DataFrame
        A pandas DataFrame with the Pb/Pb data.  Must have columns '206Pb/204Pb' and '207Pb/204Pb'
    color : string
        The color of the line.
    label : string
        The label of the line.
        
    Returns
    -------
    figure 
        A figure with a uranogenic Pb/Pb curve.

    Raises
    ------
    None.
    """
    spl = CubicSpline(df['206Pb/204Pb'], df['207Pb/204Pb'])
    y_smooth = spl(df['206Pb/204Pb'])
    fig.line(df['206Pb/204Pb'], y_smooth, color=color, legend_label=label)
    return fig

def plot_thorogenicPb_curve(fig,df,color,label):
    """Plot a thorogenic Pb/Pb curve.  206Pb/204Pb vs. 208Pb/204Pb.
    Parameters
    ----------
    fig : figure
        The figure you want to plot the diagram on.
    df : DataFrame
        A pandas DataFrame with the Pb/Pb data.  Must have columns '206Pb/204Pb' and '208Pb/204Pb'
    color : string
        The color of the line.
    label : string
        The label of the line.
        
    Returns
    -------
    figure 
        A figure with a thorogenic Pb/Pb curve.

    Raises
    ------
    None.
    """
    spl = CubicSpline(df['206Pb/204Pb'], df['208Pb/204Pb'])
    y_smooth = spl(df['206Pb/204Pb'])
    fig.line(df['206Pb/204Pb'], y_smooth, color=color, legend_label=label)
    return fig

def calc_epsilon(sample_parent, model_parent):
    """Calculates the epsilon value, given the sample parent and model parent composition.
    Parameters
    ----------
    sample_parent : float
        Measured parent composition.        
    model_parent : float
        Model parent composition. (Typically CHUR or DM for example.)
        
    Returns
    -------
    float
        The epsilon value.

    Raises
    ------
    None.
    """
    return ((sample_parent - model_parent) / model_parent)  * 10**4

def calc_delta(sample_parent, model_parent):
    """Calculates the delta value, given the sample parent and model parent composition.
    Parameters
    ----------
    sample_parent : float
        Measured parent composition.        
    model_parent : float
        Model parent composition. (Typically CHUR or DM for example.)
        
    Returns
    -------
    float
        The delta value.

    Raises
    ------
    None.
    """
    return ((sample_parent - model_parent) / model_parent)  * 10**3

def calc_model_age(t2_parent, t1_parent, t2_daughter, t1_daughter, decay_constant):
    """Calculate the model age of a sample. NOTE: THIS FUNCTION IS UNDER TESTING AND MAY NOT GIVE THE CORRECT ANSWER!
    Parameters
    ----------
    t2_parent : float
        The measured parent composition.
    t1_parent : float
        The model parent composition.
    t2_daughter : float
        The measured daughter composition.        
    t1_daughter : float
        The model daughter composition.
    decay_constant : float
        The decay constant of the isotope.
        
    Returns
    -------
    float
        The calculated model age of the sample.

    Raises
    ------
    None.
    """
    parent = t1_parent - t2_parent
    daughter = t1_daughter - t2_daughter
    return np.log((parent/daughter) + 1) / decay_constant / 1000000

