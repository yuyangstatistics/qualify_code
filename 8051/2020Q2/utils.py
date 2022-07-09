import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.model_selection import RandomizedSearchCV, GridSearchCV
import os

def connected_groups(adj):
    """Used to find the connected subgroups, given an adjacency matrix.
    The adjacency matrix can come from either correlation or shap interaction values.
    Reference: https://stackoverflow.com/questions/57825109/find-how-many-connected-groups-of-nodes-in-a-given-adjacency-matrix
    Parameters
    ----------
    adj: array, 2-dimensional 
        The adjacency matrix.
    Returns
    -------
    groups: dictionary
        A dictionary containing the subgroups.
        Key is the index of the subgroup. 
        Value is a list of indices of sensors inside the subgroup.
    
    others: list
        A list of indices of sensors which don't belong to any subgroup.
    """
    n = len(adj)
    nodes_to_check = set([i for i in range(n)])
    groups = {}
    count = 0
    others = []  # features without any grouping structue
    while nodes_to_check:
        node = nodes_to_check.pop()
        adjacent = adj[node]
        other_group_members = set()
        for i in nodes_to_check:
            if adjacent[i]:
                other_group_members.add(i)
        nodes_to_check -= other_group_members
        if other_group_members:
            count += 1
            other_group_members.add(node)
            groups[count] = list(other_group_members)
        else:
            others.append(node)

    return groups, others

def show_corr(X, cutoff = 0.9, img_dir = 'images'):
    """Investigate the correlation among the predictors.
        (1) Display correlation heatmap
        (2) Show the correlated groups in the predictors.
    
    Parameters
    ----------
    X: the predictor dataframe

    cutoff: range(0, 1), default = 0.9
        The cutoff to determine closely related or not.
    
    img_dir: str
        The directory to save the heatmap image.
    
    Returns
    -------
    No explicit return. The correlated groups will be printed out.
    """
    
    if not os.path.exists(img_dir):
        os.makedirs(img_dir)

    # plot correlation heatmap
    corr = X.corr()
    sns.set(rc={'figure.figsize':(11.7,8.27)})
    ax = sns.heatmap(
        corr, 
        vmin=-1, vmax=1, center=0,
        cmap=sns.diverging_palette(20, 220, n=200),
        square=True
    )
    ax.set_xticklabels(
        ax.get_xticklabels(),
        rotation=45,
        horizontalalignment='right'
    )
    plt.title('Heatmap of Variable Correlation')
    plt.savefig(f"{img_dir}/heatmap.pdf")
    print(f"The image has been saved to '{img_dir}/heatmap.pdf'.")
    plt.show()
    sns.set(rc={'figure.figsize':(6, 4)})
    
    ## obtain correlated groups and print them out
    corr_values = corr.abs().to_numpy()
    # get adjacency matrix
    corr_adj = (corr_values > cutoff) * 1  # rule-of-thumb is to use 0.8 as the cutoff
    # get subgroups
    subgroups, others = connected_groups(corr_adj)
    print(f"Grouping:")
    for key, value in subgroups.items():
        print(f"{key}: {list(X.columns[sorted(value)])}")
    print(f"Others: {list(X.columns[sorted(others)])}")




def search_pipeline(X_train_data, X_test_data, y_train_data, y_test_data, 
                       model, param_grid, cv=10, scoring_fit='neg_mean_squared_error',
                       do_probabilities = False, search_mode = 'GridSearchCV', n_iterations = 0):
    """ Parameter tuning pipeline for XGBoost, LightGBM, and general sklearn models.

    Reference: https://mlfromscratch.com/gridsearch-keras-sklearn/#/
    """
    fitted_model = None
    
    if(search_mode == 'GridSearchCV'):
        gs = GridSearchCV(
            estimator=model,
            param_grid=param_grid, 
            cv=cv, 
            n_jobs=-1, 
            scoring=scoring_fit,
            verbose=2
        )
        fitted_model = gs.fit(X_train_data, y_train_data)

    elif (search_mode == 'RandomizedSearchCV'):
        rs = RandomizedSearchCV(
            estimator=model,
            param_distributions=param_grid, 
            cv=cv,
            n_iter=n_iterations,
            n_jobs=-1, 
            scoring=scoring_fit,
            verbose=2
        )
        fitted_model = rs.fit(X_train_data, y_train_data)
    
    
    if(fitted_model != None):
        if do_probabilities:
            pred = fitted_model.predict_proba(X_test_data)
        else:
            pred = fitted_model.predict(X_test_data)
            
        return fitted_model, pred