import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
import seaborn as sns
import pandas as pd
from pybib.calc import *
import os
import matplotlib.font_manager as fm
from matplotlib.gridspec import GridSpec
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_squared_error, r2_score

def perform_linear_regression(x, y):

    x = np.array(x).reshape(-1, 1)  # Reshape for sklearn
    y = np.array(y)

    model = LinearRegression()
    model.fit(x, y)

    y_pred = model.predict(x)
    r2 = r2_score(y, y_pred)
    coefficients = model.coef_
    intercept = model.intercept_
    return r2, coefficients,intercept


FONT_LABEL = './font/Quicksand_Book.otf'

# Registrar a fonte no FontManager
fm.fontManager.addfont(FONT_LABEL)
prop_label = fm.FontProperties(fname=FONT_LABEL)

FONT_TICKS = './font/HostGrotesk-Regular.ttf'

# Registrar a fonte no FontManager
fm.fontManager.addfont(FONT_TICKS)
prop_ticks = fm.FontProperties(fname=FONT_TICKS)

def calculate_average_std(df, column_index,filename,value = 0):
    result = []
    for file, N in zip(df['file_name'], df['network_size']):
        data = np.loadtxt(f'{filename}/{file}').T
        if(value != 0):
            data[column_index] = np.abs(data[column_index] - value)/value
            media, std = np.mean(data[column_index]), np.std(data[column_index])
        else:
            
            media, std = np.mean(data[column_index]), np.std(data[column_index])
        result.append([N, media, std])
    return np.array(result)

def calculate_average_std_all_probabilities(df, column_index,filename):
    result = []
    for prob in df['probability'].unique():
        df_prob = df[df['probability'] == prob]
        for file, N in zip(df_prob['file_name'], df_prob['network_size']):
            data = np.loadtxt(f'{filename}/{file}').T
            media, std = np.mean(data[column_index]), np.std(data[column_index], ddof=1)
            result.append([N, prob, media, std])
    return np.array(result)

def calculate_t_test(df, column_index, y_line):
    result = []
    for file, N in zip(df['file_name'], df['network_size']):
        data = np.loadtxt(f'./output/blogs/{file}').T
        t_stat, p_value = stats.ttest_1samp(data[column_index], y_line)
        result.append([N, t_stat, p_value])
    return np.array(result)

def plot_compara_p(pvalues,df, column, label_x, label_y,filename, title="", y_line = {}, test_t=False, ispercentage=False, ax=None,linear_regression = False,value = 0):
    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 8))
    colors = ["#f1faee", "tab:cyan", "tab:brown", "tab:gray", "tab:olive", "tab:purple", "tab:brown", "tab:pink", "tab:gray", "tab:olive", "tab:cyan"]

    check = True 
    colors = np.array(colors)
    colors = colors[np.argsort(np.array(pvalues))]
    pvalues = np.sort(pvalues)
    result_prob = []
    for p,color in zip(pvalues,colors):

        df_prob = df[df['probability'] == p]
        df_prob = df_prob.sort_values(by='network_size')
        if(label_y == "Loss"):
            result_prob = calculate_average_std(df_prob, 0,filename,value)
        else:
            result_prob = calculate_average_std(df_prob, column,filename)
        if ispercentage:
            result_prob[:, 1] *= 100  # Convert to percentage
            result_prob[:, 2] *= 100  # Convert to percentage

        ax.scatter(result_prob[:, 0], result_prob[:, 1], label=f'p = {np.round(p,2)}', c=color, edgecolors='k', zorder=2)
        ax.errorbar(result_prob[:, 0], result_prob[:, 1], yerr=result_prob[:, 2], fmt='o', capsize=4, zorder=1, c='#432818')
        if ispercentage:
            ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda y, _: f'{y:.0f}%'))

        if(linear_regression):
            print(len(result_prob[:, 0]))
            y_ = np.log10(result_prob[:, 1])
            x_ = np.log10(result_prob[:, 0])
            r2,a,b = perform_linear_regression(x_[20:], y_[20:])
            print(r2,a,b)
    if len(y_line) != 0 and check:
        linestyles = [':','--', '-.']
        for key,linestyle in zip(y_line, linestyles):

            ax.axhline(y=y_line[key], color='black', linestyle=linestyle, label=key)
        check = False
        if test_t:
            resultado_test_t = calculate_t_test(df_prob, column, y_line)
            for i in resultado_test_t:
                if i[2] < 0.05:
                    print("Rejeitamos a hipótese nula. Há evidências de diferença significativa.")
                else:
                    print("Não rejeitamos a hipótese nula. Não há evidências de diferença significativa.")
    
    ax.set_xlabel(label_x, fontsize=30, fontproperties=prop_label)
    ax.set_ylabel(label_y, fontsize=30, fontproperties=prop_label)
    ax.tick_params(axis='both', which='major', labelsize=20)
    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontproperties({'size': 20, 'family': prop_ticks.get_name()})
    if title:
        ax.set_title(title)
    if(len(pvalues) > 2):
        ax.legend(loc='center right',bbox_to_anchor=(1.5, 0.7), fontsize=20, scatterpoints=1, markerscale=2,prop={'size': 20, 'family': prop_ticks.get_name()})
    else:
        ax.legend(loc='best', fontsize=20, scatterpoints=1, markerscale=2,prop={'size': 20, 'family': prop_ticks.get_name()})
    ax.grid(linestyle=':', linewidth=0.5, zorder=0)

# Call the function

def plot_compara_size(
        df, 
        column, 
        label_x, 
        label_y,
        filename, 
        title="", 
        y_line={}, 
        legend="", 
        ispercentage=False, 
        ax=None,
        has_line = True,
        linear_regression = False,
    ):
    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 8), dpi=300)
    sizes = df['network_size'].unique()

    sizes = np.sort(sizes)
    first = sizes[0]
    first_ = 10*sizes[0]
    sizes = list(sizes[sizes%first_ == 0]) + list([first])
    sizes = np.sort(sizes)
    colors = ["tab:red", "tab:blue", "tab:green", "tab:orange", "tab:purple", "tab:brown", "tab:pink", "tab:gray", "tab:olive", "tab:cyan"]
    colors = np.array(colors)
    for i,cor in zip(sizes,colors):

        df_size = df[df['network_size'] == i]
        df_size = df_size.sort_values(by='probability')
        result = calculate_average_std_all_probabilities(df_size, column,filename)

        if ispercentage:
            result[:, 2] *= 100  # Convert to percentage
        ax.scatter(result[:, 1], result[:, 2], edgecolors='k',color = cor, label=f'N* = {i}' if i < 500 else f'N = {i}', zorder=2)
        ax.errorbar(result[:, 1], result[:, 2], yerr=result[:, 3], fmt='o', capsize=4, zorder=1, c='#432818')
        if ispercentage:
            ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda y, _: f'{y:.0f}%'))
        if(df_size.shape[0] != 101):
            if(has_line):
                ax.plot(result[:, 1], result[:, 2],color = cor, zorder=1)
        if((linear_regression) & (i > 5)):
            print(i)
            y_ = np.log10(result[1:, 2])
            x_ = np.log10(result[1:, 1])
            if((df_size.shape[0] == 101)):
                r2,a,b = perform_linear_regression(x_[8:], y_[8:])
                ax.plot(result[1:, 1],result[1:, 1]**a*10**b,c ='yellow')
            else:
                r2,a,b = perform_linear_regression(x_, y_)
                ax.plot(result[1:, 1],result[1:, 1]**a*10**b,c='yellow')
            print(r2,a,b)
        
    if len(y_line) != 0 :
        linestyles = [':','--', '-.']
        for key,linestyle in zip(y_line, linestyles):

            ax.axhline(y=y_line[key], color='black', linestyle=linestyle, label=key)
    ax.set_xlabel(label_x, fontsize=30, fontproperties=prop_label)
    ax.set_ylabel(label_y, fontsize=30, fontproperties=prop_label)
    ax.tick_params(axis='both', which='major', labelsize=20)
    if title != "":
        ax.set_title(title)
    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontproperties({'size': 20, 'family': prop_ticks.get_name()})
    
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.85), scatterpoints=1, markerscale=2, prop={'size': 20, 'family': prop_ticks.get_name()})
    ax.grid(linestyle=':', linewidth=0.5, zorder=0)

def plot_test(
        df, 
        column, 
        filename, 
    ):
    sizes = df['network_size'].unique()

    sizes = np.sort(sizes)
    print(sizes)
    y = []
    for i in sizes:
        
        print(i)
        df_size = df[df['network_size'] == i]
        df_size = df_size.sort_values(by='probability')
        result = calculate_average_std_all_probabilities(df_size, column,filename)
        y_ = np.log10(result[1:, 1])
        x_ = np.log10(result[1:, 2])
        if((df_size.shape[0] == 101)):
            r2,a,b = perform_linear_regression(x_[20:], y_[20:])
            y.append(a[0])
        else:
            r2,a,b = perform_linear_regression(x_, y_)
            y.append(a[0])
    return sizes,y