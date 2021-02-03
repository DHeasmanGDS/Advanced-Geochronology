"""
Author                      : Drew Heasman
Date(last updated)          : 28 January 2021
Description : Advanced Geochronology code.
"""

import geochron_apps as geochron
import os
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

os.chdir('C:/USask Python/geochronology')
pd.options.display.max_columns = 2000
pd.set_option('display.float_format', lambda x: '%.11f' % x)

####### Concordia Diagram U-Pb Calculations
figsize = (10,12)
geochron.get_figure(figsize=figsize,
                    title="Wetherill Concordia Diagram",
                    xlabel="207Pb/235U",
                    ylabel="206Pb/238U",
                    xlim=[0,100],
                    ylim=[0,1.2])
geochron.plot_wetherill_concordia()
geochron.plot_wetherill_concordia(lambda238=1.55*10**-10,
                                        lambda235=9.72*10**-10,
                                        linecolor="g",
                                        show_ages=False,
                                        label="U Decay Constant 2020?")
plt.legend()
plt.show()
#
geochron.get_figure(figsize=figsize,
                    title="Tera-Wasserburg Condordia Diagram",
                    xlabel="238U/206Pb",
                    ylabel="207Pb/206Pb",
                    xlim=[0,75],
                    ylim=[0,0.62])

geochron.plot_tera_wasserburg_concordia()
geochron.plot_tera_wasserburg_concordia(lambda238=1.54*10**-10,
                                        lambda235=9.72*10**-10,
                                        linecolor="g",
                                        show_ages=False,
                                        label="U Decay Constant 2020?")
plt.legend()
plt.show()

####### Pb-Pb Calculations

print(geochron.calc_present_day_comp_Pb(16.0, np.array([0.1,2.0,9.74,50]), 2000))

####### Rb-Sr Calculations

print(geochron.calc_present_day_comp_Sr(0.704, 500, 2000))
print(geochron.calc_age_Sr(0.704,1000,0.73))
print(geochron.calc_age_Sr(0.0,1000,0.73))
print(geochron.calc_present_day_comp_Sr(0.704,500,500))
