# coding: utf-8

class CalculusError(Exception):
    def __init__(self,*args,**kwargs):
        Exception.__init__(self,*args,**kwargs)

class ArgumentTypeError(CalculusError):
    pass
