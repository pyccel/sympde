# coding: utf-8

class UnconsistentError(Exception):
    def __init__(self,*args,**kwargs):
        Exception.__init__(self,*args,**kwargs)

class UnconsistentLhsError(UnconsistentError):
    pass

class UnconsistentRhsError(UnconsistentError):
    pass

class UnconsistentBCError(UnconsistentError):
    pass

class UnconsistentArgumentsError(UnconsistentError):
    pass

class UnconsistentLinearExpressionError(UnconsistentError):
    pass
