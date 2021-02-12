import pandas as pd
from sklearn.base import BaseEstimator, TransformerMixin



class GetOrographyError(BaseEstimator, TransformerMixin):
    """Prepare orography error from station altitude an model grid-point altitude."""

    def fit(self, X, y=None):
        return self

    def transform(self, X, y=None):
        return pd.DataFrame({'orog_error': X.alt - X.orog})

    def get_feature_names(self):
        return ['orog_error']