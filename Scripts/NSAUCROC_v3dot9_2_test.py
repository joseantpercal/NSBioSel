### Libraries ###
import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.tree import DecisionTreeClassifier
from sklearn.metrics import fbeta_score
from sklearn.metrics import precision_score
from sklearn.metrics import confusion_matrix
from sklearn.metrics import roc_auc_score
import itertools
import time
import os

########################################################################################################

### My functions ###

def myAUCROCFeatureSelection(adata, df_dummies, cl, n_top_genes):

    start_time = time.time()  # Registrar el tiempo de inicio

    feature_scores = []

    x_train = adata.to_df()
    y_train = df_dummies[cl]

    for col in x_train.columns:       
        auc = roc_auc_score(y_train, x_train[col]) 
        if auc < 0.5: # in case the feature is negatively correlated with the target
            auc = roc_auc_score(y_train, -x_train[col]) 
        feature_scores.append((col, auc))

    # Convertimos la lista de tuplas en un diccionario para facilitar la creación de la Serie de Pandas
    feature_dict = {gene: score for gene, score in feature_scores}

    # Creamos la Serie de Pandas utilizando el diccionario y especificando el nombre del índice
    df_scores = pd.Series(feature_dict.values(), index=adata.var_names).sort_values(ascending=False)[:n_top_genes]

    end_time = time.time()  # Registrar el tiempo de finalización
    elapsed_time = end_time - start_time  # Calcular el tiempo transcurrido
    
    print("Tiempo de ejecución: {:.2f} segundos".format(elapsed_time))

    return df_scores

## construct decision tree for each gene and evaluate the fbeta score in all combinations ==> outputs markers with max fbeta, and all scores
def myDecisionTreeEvaluation(adata, df_dummies, cl, genes_eval, beta):
    dict_pred = {}
    dict_thresholds = {}  # Diccionario para almacenar los thresholds
    for i in genes_eval:
#         x_train = adata[:,i].X
        x_train = adata[:,i].to_df()
        y_train = df_dummies[cl]
        tree_clf = DecisionTreeClassifier(max_leaf_nodes=2)
        tree_clf = tree_clf.fit(x_train, y_train) 
        dict_pred[i] = tree_clf.apply(x_train)-1
        dict_thresholds[i] = tree_clf.tree_.threshold[0]  # Guardar el threshold
    df_pred = pd.DataFrame(dict_pred)
    
    combs = []
    for L in range(1, len(genes_eval)+1):
        els = [list(x) for x in itertools.combinations(genes_eval, L)]
        combs.extend(els)
    
    dict_scores = {} 
    for ii in combs:
        y_true = df_dummies[cl]
        y_pred = df_pred[ii].product(axis=1)
        fbeta = fbeta_score(y_true, y_pred, average='binary', beta=beta)
        ppv = precision_score(y_true, y_pred, average='binary', zero_division=0)
        tn, fp, fn, tp = confusion_matrix(y_true, y_pred).ravel()
        dict_scores['&'.join(ii)] = fbeta, ppv, tn, fp, fn, tp
    df_scores = pd.DataFrame(dict_scores)
        
    ## find which combination has the max fbeta
    idx_max = df_scores.idxmax(axis=1)[0] #[0] is fbeta
    markers = idx_max.split('&')
    scores = df_scores[idx_max]
    score_max = scores[0]

    # Obtener el threshold correspondiente a los marcadores seleccionados
    thresholds = [dict_thresholds[marker] for marker in markers]

    return markers, scores, score_max, thresholds

## construct decision tree for each gene and evaluate the fbeta score in all combinations ==> outputs markers with max fbeta, and all scores
def myDecisionTreeEvaluationTest(adata, adata_test, cluster_header, NSAUCROC_results, beta, output_folder, outputfilename_prefix):
    
    # Datos de ejemplo
    genes_data = NSAUCROC_results['NSAUCROC_markers']
    cluster_names = NSAUCROC_results['clusterName']

    # Crear un DataFrame vacío
    merged_data = pd.DataFrame(columns=['ClusterName', 'Genes'])

    # Iterar sobre los datos de los genes y los nombres de los clústeres
    for i, (genes, cluster_name) in enumerate(zip(genes_data, cluster_names)):
        # Agregar una nueva fila al DataFrame con el nombre del clúster y los genes
        merged_data.loc[i] = [cluster_name, genes]

    #Preprocess data
    adata.X = adata.to_df()
    adata_test.X = adata_test.to_df()
    ## categorial cluster labels
    adata.obs[cluster_header] = adata.obs[cluster_header].astype('category')
    adata_test.obs[cluster_header] = adata_test.obs[cluster_header].astype('category')
    ## dummy/indicator for one vs. all
    df_dummies = pd.get_dummies(adata.obs[cluster_header]) #cell-by-cluster
    df_dummies_test = pd.get_dummies(adata_test.obs[cluster_header]) #cell-by-cluster

    # DataFrame para almacenar los resultados
    results_df = pd.DataFrame(columns=['clusterName', 'clusterSize', 'f_score', 'PPV', 'TN', 'FP', 'FN', 'TP', 'marker_count', 'NSAUCROC_markers', 'threshold'])

    dict_pred_test = {}
    ct = 0
    for cl in merged_data.ClusterName:
        thresholds = []  # Lista para almacenar los thresholds
        for i in merged_data.Genes[ct]:
    #         x_train = adata[:,i].X
            x_train = adata[:,i].to_df()
            y_train = df_dummies[cl]
            x_test = adata_test[:,i].to_df()
            tree_clf = DecisionTreeClassifier(max_leaf_nodes=2)
            tree_clf = tree_clf.fit(x_train, y_train)
            threshold = tree_clf.tree_.threshold[0]
            thresholds.append(threshold)  # Añadir el threshold a la lista 
            dict_pred_test[i] = tree_clf.apply(x_test)-1
        df_pred_test = pd.DataFrame(dict_pred_test)
        
        gene_list = merged_data['Genes'][ct]
                
        dict_scores_test = {} 
        y_true = df_dummies_test[cl]
        y_pred = df_pred_test[gene_list].product(axis=1)
        fbeta = fbeta_score(y_true, y_pred, average='binary', beta=beta)
        ppv = precision_score(y_true, y_pred, average='binary', zero_division=0)
        tn, fp, fn, tp = confusion_matrix(y_true, y_pred).ravel()
        dict_scores_test['f_score'] = fbeta
        dict_scores_test['PPV'] = ppv
        dict_scores_test['TN'] = tn
        dict_scores_test['FP'] = fp
        dict_scores_test['FN'] = fn
        dict_scores_test['TP'] = tp
        dict_scores_test['marker_count'] = len(gene_list)
        dict_scores_test['NSAUCROC_markers'] = str(gene_list)
        
         # Crear DataFrame temporal con la fila actual
        temp_df = pd.DataFrame({'clusterName': cl, 'clusterSize': tp + fn, **dict_scores_test}, index=[0])

        # Añadir el threshold al diccionario de resultados
        temp_df['threshold'] = [thresholds]
        
        # Concatenar el DataFrame temporal con el DataFrame de resultados
        results_df = pd.concat([results_df, temp_df], ignore_index=True)

        ct += 1    

    # Calcular el promedio de f_score y PPV
    f_score_mean = results_df['f_score'].mean()
    PPV_mean = results_df['PPV'].mean()
    cluster_mean = int(results_df['clusterSize'].mean())
    TN_mean = int(results_df['TN'].mean())
    FP_mean = int(results_df['FP'].mean())
    FN_mean = int(results_df['FN'].mean())
    TP_mean = int(results_df['TP'].mean())
    
    # Construir un DataFrame para la nueva fila
    new_row = pd.DataFrame({
        'clusterName': ['Average'], 
        'clusterSize': [cluster_mean], 
        'f_score': [f_score_mean], 
        'PPV': [PPV_mean], 
        'TN': [TN_mean],
        'FP': [FP_mean],
        'FN': [FN_mean],
        'TP': [TP_mean],
        'marker_count': ['-'], 
        'NSAUCROC_markers': ['-'],
        'threshold': ['-']
    })  

    output_test = 'evaluation_test'
    # Concatenar el DataFrame de la nueva fila con el DataFrame results_df
    results_df = pd.concat([results_df, new_row], ignore_index=True)
    results_df.to_csv(f"{output_folder}{outputfilename_prefix}_{output_test}_results.csv", index=False)

    return results_df

## construct decision tree for each gene and evaluate the fbeta score in all combinations ==> outputs markers with max fbeta, and all scores
def myDecisionTreeEvaluationTestCombined(adata, adata_test, cluster_header, NSAUCROC_results, beta, output_folder, outputfilename_prefix):

    # Datos de ejemplo
    genes_data = NSAUCROC_results['NSAUCROC_markers']
    cluster_names = NSAUCROC_results['clusterName']

    # Crear un DataFrame vacío
    merged_data = pd.DataFrame(columns=['ClusterName', 'Genes'])

    # Iterar sobre los datos de los genes y los nombres de los clústeres
    for i, (genes, cluster_name) in enumerate(zip(genes_data, cluster_names)):
        # Agregar una nueva fila al DataFrame con el nombre del clúster y los genes
        merged_data.loc[i] = [cluster_name, genes]

    #Preprocess data
    adata.X = adata.to_df()
    adata_test.X = adata_test.to_df()
    ## categorial cluster labels
    adata.obs[cluster_header] = adata.obs[cluster_header].astype('category')
    adata_test.obs[cluster_header] = adata_test.obs[cluster_header].astype('category')
    ## dummy/indicator for one vs. all
    df_dummies = pd.get_dummies(adata.obs[cluster_header]) #cell-by-cluster
    df_dummies_test = pd.get_dummies(adata_test.obs[cluster_header]) #cell-by-cluster

    # DataFrame para almacenar los resultados
    results_df = pd.DataFrame(columns=['clusterName', 'clusterSize', 'f_score', 'PPV', 'TN', 'FP', 'FN', 'TP', 'marker_count', 'NSAUCROC_markers', 'threshold'])

    cluster_gene_pairs = [(row['ClusterName'], row['Genes']) for index, row in merged_data.iterrows()]
    for cl, genes in cluster_gene_pairs: 
        thresholds = []  # Lista para almacenar los thresholds
        x_train = adata[:,genes].to_df()
        y_train = df_dummies[cl]
        x_test = adata_test[:,genes].to_df()
        tree_clf = DecisionTreeClassifier(max_leaf_nodes=2)
        tree_clf = tree_clf.fit(x_train, y_train) 
        # Obtener thresholds de todos los nodos de división
        thresholds.extend(tree_clf.tree_.threshold[tree_clf.tree_.threshold != -2])  # Filtrar umbrales válidos
        y_pred = tree_clf.apply(x_test)-1
                        
        dict_scores_test = {} 
        y_true = df_dummies_test[cl]
        fbeta = fbeta_score(y_true, y_pred, average='binary', beta=beta)
        ppv = precision_score(y_true, y_pred, average='binary', zero_division=0)
        tn, fp, fn, tp = confusion_matrix(y_true, y_pred).ravel()
        dict_scores_test['f_score'] = fbeta
        dict_scores_test['PPV'] = ppv
        dict_scores_test['TN'] = tn
        dict_scores_test['FP'] = fp
        dict_scores_test['FN'] = fn
        dict_scores_test['TP'] = tp
        dict_scores_test['marker_count'] = len(genes)
        dict_scores_test['NSAUCROC_markers'] = str(genes)
        
         # Crear DataFrame temporal con la fila actual
        temp_df = pd.DataFrame({'clusterName': cl, 'clusterSize': tp + fn, **dict_scores_test}, index=[0])

        # Añadir los thresholds al DataFrame temporal
        temp_df['threshold'] = [thresholds] if thresholds else ['-']
        
        # Concatenar el DataFrame temporal con el DataFrame de resultados
        results_df = pd.concat([results_df, temp_df], ignore_index=True)

    # Calcular el promedio de f_score y PPV
    f_score_mean = results_df['f_score'].mean()
    PPV_mean = results_df['PPV'].mean()
    cluster_mean = int(results_df['clusterSize'].mean())
    TN_mean = int(results_df['TN'].mean())
    FP_mean = int(results_df['FP'].mean())
    FN_mean = int(results_df['FN'].mean())
    TP_mean = int(results_df['TP'].mean())
    
    # Construir un DataFrame para la nueva fila
    new_row = pd.DataFrame({
        'clusterName': ['Average'], 
        'clusterSize': [cluster_mean], 
        'f_score': [f_score_mean], 
        'PPV': [PPV_mean], 
        'TN': [TN_mean],
        'FP': [FP_mean],
        'FN': [FN_mean],
        'TP': [TP_mean],
        'marker_count': ['-'], 
        'NSAUCROC_markers': ['-'],
        'threshold': ['-']
    })

    output_test = 'evaluation_combined_test'
    # Concatenar el DataFrame de la nueva fila con el DataFrame results_df
    results_df = pd.concat([results_df, new_row], ignore_index=True)
    results_df.to_csv(f"{output_folder}{outputfilename_prefix}_{output_test}_results.csv", index=False)
                
    return results_df

def myDecisionTreeEvaluationTestMove(adata, adata_test, cluster_header, NSAUCROC_results, beta, output_folder, outputfilename_prefix, coef=0.5):
    
    # Datos de ejemplo
    genes_data = NSAUCROC_results['NSAUCROC_markers']
    cluster_names = NSAUCROC_results['clusterName']

    # Crear un DataFrame vacío
    merged_data = pd.DataFrame(columns=['ClusterName', 'Genes'])

    # Iterar sobre los datos de los genes y los nombres de los clústeres
    for i, (genes, cluster_name) in enumerate(zip(genes_data, cluster_names)):
        # Agregar una nueva fila al DataFrame con el nombre del clúster y los genes
        merged_data.loc[i] = [cluster_name, genes]

    #Preprocess data
    adata.X = adata.to_df()
    adata_test.X = adata_test.to_df()
    ## categorial cluster labels
    adata.obs[cluster_header] = adata.obs[cluster_header].astype('category')
    adata_test.obs[cluster_header] = adata_test.obs[cluster_header].astype('category')
    ## dummy/indicator for one vs. all Random Forest model
    df_dummies = pd.get_dummies(adata.obs[cluster_header]) #cell-by-cluster
    df_dummies_test = pd.get_dummies(adata_test.obs[cluster_header]) #cell-by-cluster

    # DataFrame para almacenar los resultados
    results_df = pd.DataFrame(columns=['clusterName', 'clusterSize', 'f_score', 'PPV', 'TN', 'FP', 'FN', 'TP', 'marker_count', 'NSAUCROC_markers', 'original_threshold', 'new_threshold'])

    dict_pred_test = {}
    ct = 0
    for cl in merged_data.ClusterName:
        original_thresholds = []  # Lista para almacenar los thresholds originales
        new_thresholds = []  # Lista para almacenar los thresholds nuevos
        for i in merged_data.Genes[ct]:
            x_train = adata[:,i].to_df()
            y_train = df_dummies[cl]
            x_test = adata_test[:,i].to_df()

            # Calcular la media y la mediana de x_train
            mean_value = np.mean(x_train[x_train > 0])
            median_value = (x_train[x_train > 0]).median()
            median_value = median_value[i]
            std_value = np.std(x_test)
            std_value = std_value[i]
            
            # Entrenar el árbol de decisión
            tree_clf = DecisionTreeClassifier(max_leaf_nodes=2)
            tree_clf = tree_clf.fit(x_train, y_train) 

            # Obtener el threshold original
            original_threshold = tree_clf.tree_.threshold[0]
            original_thresholds.append(original_threshold)

            # Calcular el nuevo threshold
            new_threshold = original_threshold + coef*std_value
            new_thresholds.append(new_threshold)

            # Asignar el nuevo threshold
            tree_clf.tree_.threshold[0] = new_threshold

            dict_pred_test[i] = tree_clf.apply(x_test)-1
        df_pred_test = pd.DataFrame(dict_pred_test)
        
        gene_list = merged_data['Genes'][ct]
                
        dict_scores_test = {} 
        y_true = df_dummies_test[cl]
        y_pred = df_pred_test[gene_list].product(axis=1)
        fbeta = fbeta_score(y_true, y_pred, average='binary', beta=beta)
        ppv = precision_score(y_true, y_pred, average='binary', zero_division=0)
        tn, fp, fn, tp = confusion_matrix(y_true, y_pred).ravel()
        dict_scores_test['f_score'] = fbeta
        dict_scores_test['PPV'] = ppv
        dict_scores_test['TN'] = tn
        dict_scores_test['FP'] = fp
        dict_scores_test['FN'] = fn
        dict_scores_test['TP'] = tp
        dict_scores_test['marker_count'] = len(gene_list)
        dict_scores_test['NSAUCROC_markers'] = str(gene_list)
        
        # Crear DataFrame temporal con la fila actual
        temp_df = pd.DataFrame({'clusterName': cl, 'clusterSize': tp + fn, **dict_scores_test}, index=[0])

        # Añadir los thresholds originales y nuevos al DataFrame temporal
        temp_df['original_threshold'] = [original_thresholds]
        temp_df['new_threshold'] = [new_thresholds]
        
        # Concatenar el DataFrame temporal con el DataFrame de resultados
        results_df = pd.concat([results_df, temp_df], ignore_index=True)

        ct += 1    

    # Calcular el promedio de f_score y PPV
    f_score_mean = results_df['f_score'].mean()
    PPV_mean = results_df['PPV'].mean()
    cluster_mean = int(results_df['clusterSize'].mean())
    TN_mean = int(results_df['TN'].mean())
    FP_mean = int(results_df['FP'].mean())
    FN_mean = int(results_df['FN'].mean())
    TP_mean = int(results_df['TP'].mean())
    
    # Construir un DataFrame para la nueva fila
    new_row = pd.DataFrame({
        'clusterName': ['Average'], 
        'clusterSize': [cluster_mean], 
        'f_score': [f_score_mean], 
        'PPV': [PPV_mean], 
        'TN': [TN_mean],
        'FP': [FP_mean],
        'FN': [FN_mean],
        'TP': [TP_mean], 
        'marker_count': ['-'], 
        'NSAUCRIC_markers': ['-'],
        'original_threshold': ['-'],
        'new_threshold': ['-']
    })

    output_test = 'evaluation_test'
    # Concatenar el DataFrame de la nueva fila con el DataFrame results_df
    results_df = pd.concat([results_df, new_row], ignore_index=True)
#    results_df.to_csv(f"{output_folder}{outputfilename_prefix}_{output_test}_results.csv", index=False)

    return results_df

###################
## Main function ##
###################

def NSAUCROC(adata, cluster_header, cluster_list=None, medians_header=None,
             beta=0.2, n_top_genes=15, n_binary_genes=10, n_genes_eval=6,
             output_folder="", outputfilename_prefix = ""):
    
    ## set up outpur folder
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    
    print("Preparing data...")
    start_global_time = time.time()
    ## densify X from sparse matrix format
    adata.X = adata.to_df()
    ## categorial cluster labels
    adata.obs[cluster_header] = adata.obs[cluster_header].astype('category')
    ## dummy/indicator for one vs. all
    df_dummies = pd.get_dummies(adata.obs[cluster_header]) #cell-by-cluster
    ## get number of cluster
    n_total_clusters = len(df_dummies.columns)

    if medians_header == None:
        print("Calculating medians...")
        start_time = time.time()
        ## get dataframes for X and cluster in a column
        df_X = pd.DataFrame(adata.X, index=adata.obs_names, columns=adata.var_names) #cell-by-gene
        clusters = adata.obs[cluster_header]
        df_X_clusters = pd.concat([df_X, clusters], axis=1)
        ## get cluster medians
        cluster_medians = df_X_clusters.groupby([cluster_header]).median() #cluster-by-gene
        ## delete to free up memories
        del df_X, clusters, df_X_clusters
        print("--- %s seconds ---" % (time.time() - start_time))
    else:
        print("Getting pre-calculated medians...")
        start_time = time.time()
        cluster_medians = adata.varm[medians_header].transpose() #cluster-by-gene
        print("--- %s seconds ---" % (time.time() - start_time))
    
    ### START iterations ###
    if cluster_list == None:
        cluster_list = df_dummies.columns
    n_clusters = len(cluster_list)
    
    print ("Number of clusters to evaluate: " + str(n_clusters))
    ct = 0
    df_supp = df_markers = df_results = pd.DataFrame()
    start_time = time.time()
    for cl in cluster_list:
        ct+=1
        print(str(ct) + " out of " + str(n_clusters) + ":")

        ## cluster in iteration
        print("\t" + str(cl))
        
        ##=== reset parameters for this iteration!!! (for taking care of special cases) ===##
        n_binary_genes_cl = n_binary_genes
        n_genes_eval_cl = n_genes_eval

        ## Feature selection step: get top genes ranked by correlation between each feature and target
        top_ROC_genes = myAUCROCFeatureSelection(adata, df_dummies, cl, n_top_genes)

        ## filter out negative genes by thresholding median>0 ==> to prevent deviding by 0 in binary score calculation
        top_gene_medians = cluster_medians.loc[cl,top_ROC_genes.index]
        top_ROC_genes_positive = top_gene_medians[top_gene_medians>0]
        n_positive_genes = sum(top_gene_medians>0)
    
        ##=== special cases: ===##
        if n_positive_genes == 0:
            print("\t" + "No positive genes for evaluation. Skipped. Optionally, consider increasing n_top_genes.")
            continue

        if n_positive_genes < n_binary_genes:
            print("\t" + f"Only {n_positive_genes} out of {n_top_genes} top selected features with median > 0 will be further evaluated.")
            n_binary_genes_cl = n_positive_genes
            n_genes_eval_cl = min(n_positive_genes, n_genes_eval)
        ##===##
            
        ## Binary scoring step: calculate binary scores for all positive top genes
        binary_scores = [sum(np.maximum(0,1-cluster_medians[i]/cluster_medians.loc[cl,i]))/(n_total_clusters-1) for i in top_ROC_genes_positive.index]
        top_binary_genes = pd.Series(binary_scores, index=top_ROC_genes_positive.index).sort_values(ascending=False)

        ## Evaluation step: calculate F-beta score for gene combinations
        genes_eval = top_binary_genes.index[:n_genes_eval_cl].to_list()
        markers, scores, score_max, thresholds = myDecisionTreeEvaluation(adata, df_dummies, cl, genes_eval, beta)
        print("\t" + str(markers))
        print("\t" + "fbeta: " + str(score_max))

        ## return supplementary table as csv
        binary_genes_list = top_binary_genes.index[:n_binary_genes_cl].to_list()
        df_supp_cl = pd.DataFrame({'clusterName': cl,
                                   'binary_genes': binary_genes_list,
                                   'AUCROC_feat_imp': top_ROC_genes[binary_genes_list],
                                   'cluster_median': top_gene_medians[binary_genes_list],
                                   'binary_score': top_binary_genes[binary_genes_list]}).sort_values('binary_score', ascending=False)
        df_supp = pd.concat([df_supp,df_supp_cl]).reset_index(drop=True)
        df_supp.to_csv(f"{output_folder}{outputfilename_prefix}_supplementary.csv", index=False)

        ## return markers table as csv
        df_markers_cl = pd.DataFrame({'clusterName': cl, 'markerGene': markers, 'score': scores[0], 'thresholds': thresholds})
        df_markers = pd.concat([df_markers, df_markers_cl]).reset_index(drop=True)
        df_markers.to_csv(f"{output_folder}{outputfilename_prefix}_markers.csv", index=False)

        ## return final results as dataframe
        dict_results_cl = {'clusterName': cl,
                           'clusterSize': int(scores[4]+scores[5]),
                           'f_score': scores[0],
                           'PPV': scores[1],
                           'TN': int(scores[2]),
                           'FP': int(scores[3]),
                           'FN': int(scores[4]),
                           'TP': int(scores[5]),
                           'marker_count': len(markers),
                           'NSAUCROC_markers': [markers],
                           'thresholds': [thresholds],  # Añadir los thresholds al DataFrame de resultados
                           'binary_genes': [df_supp_cl['binary_genes'].to_list()] #for this order is the same as the supp order
                           }
        
        df_results_cl = pd.DataFrame(dict_results_cl)
        df_results = pd.concat([df_results,df_results_cl]).reset_index(drop=True)
        df_results.to_csv(f"{output_folder}{outputfilename_prefix}_results.csv", index=False)

    # Calcular el promedio de f_score y PPV
    f_score_mean = df_results['f_score'].mean()
    PPV_mean = df_results['PPV'].mean()
    cluster_mean = int(df_results['clusterSize'].mean())
    TN_mean = int(df_results['TN'].mean())
    FP_mean = int(df_results['FP'].mean())
    FN_mean = int(df_results['FN'].mean())
    TP_mean = int(df_results['TP'].mean())

    # Crear un DataFrame para la fila de promedios
    average_row = pd.DataFrame({
        'clusterName': ['Average'],
        'clusterSize': [cluster_mean],
        'f_score': [f_score_mean],
        'PPV': [PPV_mean],
        'TN': [TN_mean],
        'FP': [FP_mean],
        'FN': [FN_mean],
        'TP': [TP_mean],
        'marker_count': ['-'],
        'NSAUCROC_markers': ['-'],
        'binary_genes': ['-'],
        'thresholds': ['-']
    })

    # Concatenar la fila de promedios al DataFrame de resultados
    df_results = pd.concat([df_results, average_row], ignore_index=True)

    print(f"Saving supplementary table as...\n{output_folder}{outputfilename_prefix}_supplementary.csv")
    print(f"Saving markers table as...\n{output_folder}{outputfilename_prefix}_markers.csv")
    print(f"Saving results table as...\n{output_folder}{outputfilename_prefix}_results.csv")
    df_results.to_csv(f"{output_folder}{outputfilename_prefix}_results.csv", index=False)
    print(f"Saving final results table as...\n{output_folder}{outputfilename_prefix}_results.csv")

    end_total_time = time.time()  # Registrar el tiempo de finalización
    elapsed_total_time = end_total_time - start_global_time  # Calcular el tiempo transcurrido
    
    print("Tiempo total de la ejecución: {:.2f} segundos".format(elapsed_total_time))    ### END iterations ###    ### END iterations ###
    
    return(df_results)

########################################################################################################
