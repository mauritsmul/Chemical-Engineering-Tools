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

    def calculate_relative_volatility(antoine, T):
        """Calculates the relative volatility of the binary mixture
        
        Keyword arguments:
            antoine: list containing two lists : containing two sets of the antoine coefficients of the binary mixture
            T: float : Temperature in Kelvin
        
        Returns:
            alpha: float : Relative volatility of the binary mixture"""
        
        # Create empty list for vapor pressures
        P_vap = [10 **(i[0] - i[1] / (T + i[2])) for i in antoine]
    
        # Calculate the relative volatility
        alpha = max(P_vap)/min(P_vap)
        return alpha

    def calculate_minimum_reflux_ratio(xd, xf, alpha):
        """Calculates the minimum reflux ratio
        
        Keyword arguments:
            xd: float : Distillate molar fraction of the most volatile component
            xf: float : Feed composition of the most volatile component
            alpha: float : Relative volatility of the binary mixture
            
        Returns:
            R_min: float : Minimum reflux ratio"""
        
        # Calculate the minimum reflux ratio
        R_min = (xd/xf-alpha*((1-xd)/(1-xf)))/(alpha-1)
        return R_min

    def plot_mccabe_thiele(T, xd, xb, xf, alpha, q, reflux_ratio):
        """Plots the McCabe-Thiele diagram
        
        Keyword arguments:
            T: float : Temperature at which the distillation is done in Kelvin
            xd: float : Distillate molar fraction of the most volatile component
            xb: float : Bottom molar fraction of the most volatile component
            xf: float : Feed composition of the most volatile component
            alpha: float : Relative volatility of the binary mixture
            q: float : Quality of the feed
            reflux_ratio: float : Reflux ratio

        Returns:
            Plot of the McCabe-Thiele diagram for the given binary mixture with given input parameters

        """

        # Create xy line and equilibrium line
        x = np.linspace(0,1,100)
        equilibrium_line = alpha*x/(1+(alpha-1)*x)
        plt.plot(x, x, color='black')
        plt.plot(x, equilibrium_line, color='black')
        plt.xlabel('x')
        plt.ylabel('y')
        plt.title('McCabe-Thiele diagram')
        plt.axis([0, 1, 0, 1])
        plt.show()

        # Plot the q line
        
        # Find the intersection between the q line and the rectifying line
        
        
        # Plot the rectifying line
        
        # Plot the stripping line
        
        # Draw the stages

    def show_design_summary():
        pass
         


# Plot the McCabe-Thiele diagram
def plot_mccabe_thiele(T, xb, xf, xd, antoine, q, reflux_factor):
    
    """
    Plots the McCabe-Thiele diagram using the user input parameters

    Keyword arguments:
            T: float : Temperature at which the distillation is done in Kelvin
            xd: float : Distillate molar fraction of the most volatile component
            xb: float : Bottom molar fraction of the most volatile component
            xf: float : Feed composition of the most volatile component
            alpha: float : Relative volatility of the binary mixture
            q: float : Quality of the feed
            reflux_ratio: float : Reflux ratio

    Returns:
            Plot of the McCabe-Thiele diagram for the given binary mixture with given input parameters
    """

    alpha = Ideal_Distillation.calculate_relative_volatility(antoine, T)
    reflux_ratio = reflux_factor*Ideal_Distillation.calculate_minimum_reflux_ratio(xd, xf, alpha)
    Ideal_Distillation.plot_mccabe_thiele(T, xd, xb, xf, alpha, q, reflux_ratio)


"""INPUT PARAMETERS FROM USER"""
# Input parameters
T = 300
xb = 0.1
xf = 0.5
xd = 0.90
antoine = [[4.16272, 1371.583, -58.496],[6.07936, 2692.187, -17.94]]
q = 1
reflux_factor = 1.5 

# Run the script
plot_mccabe_thiele(T, xb, xf, xd, antoine, q, reflux_factor)
