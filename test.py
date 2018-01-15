# -*- coding: utf-8 -*-
"""
Created on Thu Jan 11 17:14:56 2018

@author: Larry
"""

class p():
    def __init__(self, p):
        self.p = p
        
    def __str__(self):
        return 'pwep'
    def __cmp(self):
        return 2
    def __float__(self):
        return 2.2
    def __add__(self, b):
        return self.p + 3
    
a = p(2)

print (a)
print (a + 2)