"""
Author                      : Drew Heasman
Date(last updated)          : 30 January 2021
Description : Advanced Geochronology code.
"""

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

def get_figure(xlim, ylim, figsize=(3,5), xlabel="", ylabel="", title=""):
    """ Generates a blank figure for Concordia Diagrams
        
        Parameters:
           xlim (list): x range for the diagram
           ylim (list): y range for the diagram
           figsize (tuple): figure size
           xlabel (string): x axis label
           ylabel (string): y axis label
           title (string): Title for the plot
                  
        Raises:
                   
        Returns:
            No return
    """
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(1, 1, 1)
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    plt.xlabel(xlabel, fontsize=16)
    plt.ylabel(ylabel, fontsize=16)
    plt.title(title, fontsize=22)

def get_ratio_data(lambda238, lambda235, age=np.append([0.001], [np.arange(1,4568)])):
    """ Generates a DataFrame with various ratios for concordia plots. 
        
        Parameters:
           lambda238 (float): 238U decay constant, default is 1975 IUGS constant
           lambda235 (float): 235U decay constant, default is 1975 IUGS constant
                  
        Raises:
                   
        Returns:
            pandas DataFrame
    """
        
    age = age
    df = pd.DataFrame()
    df["age"] = age
    df["ratio_238U_206Pb"] = 1/(np.exp(lambda238*df.age*1000000)-1)
    df["ratio_207Pb_235U"] = np.exp(lambda235*df.age*1000000)-1
    df["ratio_206Pb_238U"] = np.exp(lambda238*df.age*1000000)-1
    df["ratio_207Pb_206Pb"] = (1/137.88)*(df.ratio_207Pb_235U/df.ratio_206Pb_238U)    
    return df

def get_ages(lambda238, lambda235, kind):
    """ Generates a plot of ages on a concordia diagram. 
        
        Parameters:
           lambda238 (float): 238U decay constant, default is 1975 IUGS constant
           lambda235 (float): 235U decay constant, default is 1975 IUGS constant
                  
        Raises:
                   
        Returns:
            No return
    """
    if kind == 'wetherill':
        x = "ratio_207Pb_235U"
        y = "ratio_206Pb_238U"
        xtext = 20
        ytext = 0
    elif kind == 'tw':
        x = "ratio_238U_206Pb"
        y = "ratio_207Pb_206Pb"
        xtext = 0
        ytext = 20
    df = get_ratio_data(lambda238, lambda235)
    ages = df.loc[
             (df['age'] == 100) |
             (df['age'] == 200) |
             (df['age'] == 500) |
             (df['age'] == 1000) |
             (df['age'] == 2000) |
             (df['age'] == 3000) |
             (df['age'] == 4000) |
             (df['age'] == 4500)]
    plt.scatter(ages[x], ages[y], marker="o", s=50, zorder=2)
    for index, row in ages.iterrows():
            label = row.age
            label = str(label) + " Ma"
            plt.annotate(label,
                (row[x], row[y]),
                 textcoords="offset points",
                 xytext=(xtext,ytext),
                 arrowprops = dict(facecolor = 'blue',
                                   arrowstyle = '->'),
                                   
                 ha='left')
    

def plot_wetherill_concordia(lambda238=1.55125*10**-10,
                   lambda235=9.8485*10**-10,
                   show_ages=True,
                   linecolor='r',
                   label="U Decay Rate 1975"):
    """ Generates a Wetherill Concodria diagram. 
        
        Parameters:
           lambda238 (float): 238U decay constant, default is 1975 IUGS constant
           lambda235 (float): 235U decay constant, default is 1975 IUGS constant
           show_ages (boolean): If True, will plot ages on the diagram
           linecolor (string): Uses matplotlib library to plot line color
                  
        Raises:
                   
        Returns:
            No return
    """
    lambda238 = lambda238
    lambda235 = lambda235

    df = get_ratio_data(lambda238, lambda235)    

    plt.plot(df["ratio_207Pb_235U"], df["ratio_206Pb_238U"],  
                 linecolor, label=label, zorder=1)
    if show_ages == True:
        get_ages(lambda238, lambda235, 'wetherill')
    
        
def plot_tera_wasserburg_concordia(lambda238=1.55125*10**-10,
                   lambda235=9.8485*10**-10,
                   show_ages=True,
                   linecolor="r",
                   label="U Decay Constant 1975"):
    """ Generates a Tera-Wasserburg Concodria diagram. 
        
        Parameters:
           lambda238 (float): 238U decay constant, default is 1975 IUGS constant
           lambda235 (float): 235U decay constant, default is 1975 IUGS constant
           show_ages (boolean): If True, will plot ages on the diagram
           linecolor (string): Uses matplotlib library to plot line color
                  
        Raises:
                   
        Returns:
            No return
    """
    df = get_ratio_data(lambda238, lambda235)    

    if show_ages == True:
        get_ages(lambda238, lambda235, 'tw')
    plt.plot(df["ratio_238U_206Pb"], df["ratio_207Pb_206Pb"],
             linecolor, label=label, zorder=1)

def calc_present_day_comp_Pb(initial,
                              mu,
                              age,
                              lambda238=1.55125*10**-10):
    """ Gives the present day composition of 206Pb/204Pb. Given a starting
        composition, 238U/204Pb(mu) values, and when minerals formed in Ma.
        
        Parameters:
            initial (float): Starting composition of mineral 206Pb/204Pb
            mu (float): 238U/204Pb ratios (often uses the mu symbol)
            age (float): When the minerals formed in Ma.
            lambda238 (float): 238U decay constant, default is 1975 IUGS constant
           
        Raises:
                   
        Returns:
            Returns the calculated 206Pb/204Pb composition in present day.
    """
    return initial + (mu * (np.exp(lambda238*age*1000000)-1))

def calc_present_day_comp_Sr(initial,
                              mu,
                              age,
                              lambda87=1.42*10**-11):
    """ Gives the present day composition of 87Sr/86Sr. Given a starting
        composition, 87Sr/86Sr(mu) values, and when minerals formed in Ma.
        
        Parameters:
            initial (float): Starting composition of mineral 87Sr/86Sr
            mu (float): 87Sr/86Sr ratios (often uses the mu symbol)
            age (float): When the minerals formed in Ma.
            lambda238 (float): Rb decay constant, default is 1975 IUGS constant
           
        Raises:
                   
        Returns:
            Returns the calculated 87Sr/87Sr composition in present day.
    """
    return initial + (mu * (np.exp(lambda87*age*1000000)-1))

def calc_age_Sr(est_initial,
                mu,
                present_day,
                lambda87=1.42*10**-11):
    """ Gives an age estimate given 87Sr/86Sr. Given a starting
        composition, 87Sr/86Sr(mu) values, and an estimate of composition at formation.
        
        Parameters:
            est_initial (float): Starting composition of mineral 87Sr/86Sr
            mu (float): 87Sr/86Sr ratios (often uses the mu symbol)
            age (float): When the minerals formed in Ma.
            lambda238 (float): Rb decay constant, default is 1975 IUGS constant
           
        Raises:
                   
        Returns:
            Returns the calculated 206Pb/204Pb composition in present day in Ma.
    """
    return np.log((present_day/mu - est_initial/mu + 1)) / lambda87 / 1000000