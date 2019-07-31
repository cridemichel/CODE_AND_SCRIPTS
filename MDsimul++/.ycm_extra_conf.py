def Settings( **kwargs ):
    #language = kwargs[ 'language' ]
    #if language=='cfamily':
    #if language == 'python':
    #MACROS=''   
    #MACLST=MACROS.split(' ')
    flags = [ '-x', 'c++', '-std=c++17', '-Wall', '-Wextra', '-I ../mdlib']
    #flags=flags+MACLST
    return {'flags': flags}
