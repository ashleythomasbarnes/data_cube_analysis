import numpy as np
import math

def mag_error(fit, fit_err, distance, distance_err):
    A, B = fit
    A_err, B_err = fit_err

    A_sqr = A **2
    B_sqr = B **2
    A_sqr_err = 2 * A * A_err
    B_sqr_err = 2 * B * B_err

    A_sqr_B_sqr = A_sqr + B_sqr
    A_sqr_B_sqr_err = (A_sqr_err **2 +  B_sqr_err**2) **0.5

    A_sqr_B_sqr_sqrt = A_sqr_B_sqr **0.5
    A_sqr_B_sqr_sqrt_err = 0.5 * A_sqr_B_sqr **-0.5 * A_sqr_B_sqr_err

    function = A_sqr_B_sqr_sqrt / distance
    function_err = function * ( (A_sqr_B_sqr_sqrt_err / A_sqr_B_sqr_sqrt) ** 2 + (distance_err / distance) **2 ) **0.5

    return function, function_err

def ang_error(fit, fit_err):
    # From https://answers.yahoo.com/question/index?qid=20110327153022AANBHtU
    A, B = fit
    A_err, B_err = fit_err

    function = math.atan(A / B)

    # A_B = (A / B)
    # A_B_err = A_B * ( ( ((A_err / A) ** 2) + ((B_err / B) **2) ) **0.5 )
    # function_err = A_B_err**2 / (1 + A_B **2) #NO!

    function_err = ( ( (B **2 * A_err **2) + (A **2 * B_err**2) ) **0.5) / (A **2 + B **2)

    return function, function_err
