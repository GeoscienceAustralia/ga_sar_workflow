# Note: Due to the fact the minimum version we support is 20191203 - we just import the proxy unmodified (as the machine generated proxy is for this version already)
import insar.gamma.generated.py_gamma_proxy
import inspect

GammaProxyBase = insar.gamma.generated.py_gamma_proxy.GammaProxy

class GammaProxy(GammaProxyBase):
    pass
