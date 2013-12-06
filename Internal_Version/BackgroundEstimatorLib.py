#$Id: BackgroundEstimatorLib.py,v 1.3 2011/10/06 10:41:50 vlasisva Exp $
def generate(env, **kw):
    if not kw.get('depsOnly',0):
        env.Tool('addLibrary', library = ['BackgroundEstimator'], package = 'BackgroundEstimator')
        pass
    env.Tool('astroLib')
    env.Tool('st_facilitiesLib')
    #env.Tool('LikelihoodLib')
    env.Tool('irfInterfaceLib')
    env.Tool('astroLib')
    env.Tool('irfLoaderLib')
    env.Tool('rootIrfLoaderLib')
    env.Tool('addLibrary', library=env['clhepLibs'])
    env.Tool('addLibrary', library=env['rootLibs'])
    #env.Tool('addLibrary', library=env['STforGRBLibs'])
    pass

def exists(env):
    return 1
