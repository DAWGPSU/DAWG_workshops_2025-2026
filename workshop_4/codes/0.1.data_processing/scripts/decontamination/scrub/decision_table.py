import sys
import numpy as np
import pandas as pd

path_before = sys.argv[1]
path_after = sys.argv[2]
path_decision = sys.argv[3]

def main(path_before, path_after, path_decision):
    otu_decontaminated = pd.read_csv(path_after, index_col=0)
    otu_before = pd.read_csv(path_before, index_col=0)
    otu_before = otu_before.loc[otu_decontaminated.index.values]

    # means - more difference, more value
    coefs = pd.DataFrame({taxa:(-otu_decontaminated[taxa].values+otu_before[taxa].values)/otu_before[taxa].values
                          for taxa in otu_decontaminated.columns.values})
    means = np.nanmean(coefs.values, axis=0)
    nans = coefs.isna().sum(axis=0)/coefs.shape[0]
    decision_table = pd.DataFrame({'Scrub means': means, 'Scrub nans':nans}, index=coefs.columns.values)
    decision_table.to_csv(path_decision)

if __name__ == '__main__':
    main(path_before, path_after, path_decision)