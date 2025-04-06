import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as sp

class Ideal_Distillation:


    """
    The Ideal_Distillation class is used to plot the McCabe-Thiele method 
    for the design of a distillation column containing an ideal binary mixture.
    
    Parameters:
    
    T: float : Temperature at which the distillation is done in Kelvin
    xd: float : Distillate molar fraction of the most volatile component
    xb: float : Bottom molar fraction of the most volatile component
    xf: float : Feed composition of the most volatile component
    antoine: [[a, b, c], [a, b, c]] : Antoine coefficients of the binary mixture 
    q: float : Quality of the feed
    R: float : Reflux factor with which the minimum reflux ratio is multiplied (Conventionally R_fac = 1.2-1.5)

    """


    def __init__(self, T: float, xd: float, xb:float, xf:float, antoine:list, q:float, R_fac:float):
        self.T = T
        self.xd = xd
        self.xb = xb
        self.xf = xf
        self.antoine = antoine
        self.q = q
        self.R_fac = R_fac

    def plot_mccabe_thiele(antoine, T, xd, xb, xf, q, reflux_factor):
        """Plots the McCabe-Thiele diagram
        
        Keyword arguments:
            T: float : Temperature at which the distillation is done in Kelvin
            xd: float : Distillate molar fraction of the most volatile component
            xb: float : Bottom molar fraction of the most volatile component
            xf: float : Feed composition of the most volatile component
            alpha: float : Relative volatility of the binary mixture
            q: float : Quality of the feed
            reflux_factor: float : Factor which is multiplied with the minimum reflux ratio

        Returns:
            Plot of the McCabe-Thiele diagram for the given binary mixture with given input parameters
            number of stages: int : Number of stages
            reflux_ratio: float : Reflux ratio
        """

        # Calculate the vapor pressures
        P_vap = [10 **(i[0] - i[1] / (T + i[2])) for i in antoine]
    
        # Calculate the relative volatility
        alpha = max(P_vap)/min(P_vap)

        # Create xy line and equilibrium line
        x = np.linspace(0,1,100)
        y = np.linspace(0,1,100)
        equilibrium_line = alpha*x/(1+(alpha-1)*x)

        # Plotting the xy line and the equilibrium line
        plt.plot(x, y, color='black')
        plt.plot(x, equilibrium_line, color='black')
        plt.xlabel('x')
        plt.ylabel('y')
        plt.title('McCabe-Thiele diagram')
        plt.axis([0, 1, 0, 1])

        # Plot the q line 
        if q == 1:
            x_values_q_line = xf*np.ones(len(y)) 
            plt.plot(x_values_q_line, y, label='q line')
        else:
            q_line = -q/(1-q)*x + xf/(1-q)
            plt.plot(x, q_line, label='q line')
        
        # Calculate the intersection point between the q line and the equilibrium line
        if q == 1: 
            intersection_q_and_VLE = xf
        else:
            intersection_q_and_VLE = sp.fsolve(lambda x: alpha*x/(1+(alpha-1)*x) - (-q/(1-q)*x + xf/(1-q)), xf)

        # Check if xd is above the intersection between the q_line and the equilibrium line
        intersection_y = alpha*intersection_q_and_VLE/(1+(alpha-1)*intersection_q_and_VLE)
        if intersection_y > xd:
            raise ValueError('Infeasible "xd" parameter, the vapor fraction cannot decrease while the liquid molar fraction is increasing')
        elif intersection_q_and_VLE < xb:
            raise ValueError('Infeasible "xb" parameter, the vapor fraction cannot increase while the liquid molar fraction is decreasing')

        # Calculate the maximum reflux ratio
        minimum_reflux = (xd/intersection_q_and_VLE-alpha*((1-xd)/(1-intersection_q_and_VLE)))/(alpha-1)

        # Set the reflux ratio
        reflux_ratio = minimum_reflux*reflux_factor
        
        # Plot the rectifying line
        rectifying_line = (reflux_ratio / (1+reflux_ratio)) * x + ( 1/(1 + reflux_ratio)) * xd
        plt.plot(x, rectifying_line, label='Rectifying')

        # Plot the stripping line
        if q == 1: 
            x_upper = xf
        else:
            x_upper = sp.fsolve(lambda x: (reflux_ratio/(1+reflux_ratio))*x+(1/(1+reflux_ratio))*xd - (-q/(1-q)*x + xf/(1-q)), xf)

        y_upper = (reflux_ratio/(1+reflux_ratio))*x_upper+(1/(1+reflux_ratio))*xd
        stripping_slope = (y_upper-xb)/(x_upper-xb)
        b = xb - stripping_slope*xb
        y_stripping = stripping_slope*x + b
        plt.plot(x, y_stripping, label='Stripping')

        # Drawing the stages
        # First we set the conditions for when the seesaw is in the stripping or rectifying section
        # Condition is True when the drawing algorithm is in the corresponding section
        # The algorithm starts in the stripping section:
        stripping_condition, rectifying_condition = True, False
        
        # Make an array in which the points can be stored
        seesaw_points_x = [np.array([xb])]
        seesaw_points_y = [np.array([xb])]

        # Set initial starting point
        x_start, y_start = np.array([xb]), np.array([xb])

        # Draw the seesaw between the operating lines and the equilibrium line
        while stripping_condition:
            # Save the VLE point
            seesaw_points_x.append(x_start)
            vle_point = alpha*x_start/(1+(alpha-1)*x_start)
            seesaw_points_y.append(vle_point)

            # Check VLE conditional
            if vle_point > y_upper:
                stripping_condition = False
                rectifying_condition = True
                x_start = sp.fsolve(lambda x: (reflux_ratio/(1+reflux_ratio))*x+(1/(1+reflux_ratio))*xd - vle_point, x0=xf)
                y_start = vle_point
                seesaw_points_x.append(x_start)
                seesaw_points_y.append(y_start)
                break 

            # Calculate the new starting point on the stripping section
            x_start = (vle_point-b)/stripping_slope
            y_start = stripping_slope*x_start + b
            seesaw_points_x.append(x_start)
            seesaw_points_y.append(y_start)

        while rectifying_condition:
            # Save the VLE point
            seesaw_points_x.append(x_start)
            vle_point = alpha*x_start/(1+(alpha-1)*x_start)
            seesaw_points_y.append(vle_point)

            # Check VLE conditional
            if vle_point > xd:
                rectifying_condition = False
                seesaw_points_y.pop(-1)
                seesaw_points_y.append(np.array([xd]))
                # Add last point
                seesaw_points_y.append(np.array([xd]))
                seesaw_points_x.append(np.array([xd]))
                break 
            
            # Calculate the new starting point on the rectifying section
            a = reflux_ratio/(1+reflux_ratio)
            b = (1/(1+reflux_ratio))*xd
            x_start = sp.fsolve(lambda x: a*x+b - vle_point, vle_point)
            seesaw_points_x.append(x_start)
            seesaw_points_y.append(vle_point)

        # Calculate the number of stages
        number_of_stages = (len(seesaw_points_x)-1)/2 - 1

        # Plot the stages calculated
        plt.plot(seesaw_points_x, seesaw_points_y)
        plt.show()
        return number_of_stages, reflux_ratio

    def show_design_summary(reflux_ratio, T, xb, xf, xd, number_of_stages):
        
        """
        Prints a summary of the distillation design to the console.

        Parameters
        ----------
        reflux_ratio : float
            The reflux ratio of the distillation.
        T : float
            The temperature at which the distillation is done in Kelvin.
        xb : float
            The bottoms composition of the most volatile component.
        xf : float
            The feed composition of the most volatile component.
        xd : float
            The distillate composition of the most volatile component.
        number_of_stages : int
            The number of stages in the distillation column.
            
        Notes
        -----
        The summary includes the feed composition, bottoms composition, distillate composition, reflux ratio, temperature, and number of stages in the distillation column.
        """
        # Make a float from the reflux_ratio array
        reflux_ratio = reflux_ratio[0]

        summary = f"""
        ================================================
                     Distillation Design Summary
        ================================================
        | Design Variable           | Value           |
        ------------------------------------------------
        | Feed Composition (xf)     | {xf:.2f}            |
        | Bottoms Composition (xb)  | {xb:.2f}            |
        | Distillate Composition(xd)| {xd:.2f}            |
        | Reflux Ratio              | {reflux_ratio:.2f}            |
        | Temperature (°C)          | {T-273:.2f} °C        |
        | Number of Trays           | {number_of_stages}             |
        ================================================
        """
        print(summary)
         

"""INPUT PARAMETERS FROM USER"""

# Input parameters

"""
T               :    float     : Temperature at which the distillation is done in Kelvin
xb              :    float     : Bottom molar fraction of the most volatile component
xf              :    float     : Feed composition of the most volatile component
xd              :    float     : Distillate molar fraction of the most volatile component
antoine         :    list      : Antoine coeffients of the binary mixture
q               :    float     : Quality of the feed
reflux_factor   :    float     : factor with which the reflux is mulitplied (conventionally ~1.2-1.5 for real world applications)
"""

T = 350
xb = 0.1
xf = 0.5
xd = 0.99
antoine = [[4.35576,	1175.581,	-2.071],[4.02832,	1268.636,	-56.199]]
q = 0.5
reflux_factor = 1.1

# Run the script
number_of_stages, reflux_ratio = Ideal_Distillation.plot_mccabe_thiele(antoine, T, xd, xb, xf, q, reflux_factor)
Ideal_Distillation.show_design_summary(reflux_ratio, T, xb, xf, xd, number_of_stages)
