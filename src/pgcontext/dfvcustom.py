from dna_features_viewer import GraphicFeature, GraphicRecord


class CustomGraphicFeature(GraphicFeature):
    def __init__(self, *args, height=1, **kwargs):
        super().__init__(*args, **kwargs)
        self.height = height

    def set_height(self, ax):
        # Custom function to set the height of the feature
        patch = super().as_patch()
        transform = ax.transData
        scaled_height = transform.transform((0, self.height))[1] - transform.transform((0, 0))[1]
        patch.set_height(scaled_height)
        return patch