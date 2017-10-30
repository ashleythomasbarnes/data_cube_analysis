import numpy as np
import math

def velo_gradients(X, a, b):
    ave_velo, ra, dec = X
    v_lsr = ave_velo + ( a * ra ) + ( b * dec )
    return v_lsr

def velo_gradients_number(a, b, distance):
    delta_v_lsr = (a**2 + b**2)**0.5 / distance
    theta_delta_v_lsr = math.atan(a / b)
    return delta_v_lsr, theta_delta_v_lsr
