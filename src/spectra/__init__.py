__all__ = ['SpectrumAnnotator', 'Booster', 'Butterfly']

def __getattr__(name):
    if name == 'SpectrumAnnotator':
        from .annotation import SpectrumAnnotator
        return SpectrumAnnotator
    elif name == 'Booster':
        from .msbooster import Booster
        return Booster
    elif name == 'Butterfly':
        from .butterfly import Butterfly
        return Butterfly
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")

def __dir__():
    return __all__

