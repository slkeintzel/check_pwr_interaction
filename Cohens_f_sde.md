# Cohens f explanation for Power Simulation for interaction effect in 2 x 2 between subjects ANOVA Design 
# Authors: Salomé Li Keintzel, Frauke Langhein and Paula G. Rosendahl (2023)

For mean structures of interaction designs the respective (residual) sd was calculated using Cohens f = 0.4 for a balanced designs.   

Following Cohen, 1977, p.371:
For our 2x2 between subjects design with Factors A and a = 2 levels and B with b = 2 levels

Cohens_f = sigma_x/sigma_e and respectively sigma_e = sigma_x/f  
with  
sigma_e being the (residual) standard deviation per cell (assumed to be equal across cells) and  
sigma_x being sqrt(sum((mean_ab - mean_a - mean_b + mean_total)^2) / (a*b)) and mean_total being the unweighed mean of cell means  

sigma_x directly translates to effect coded b3 for our three interaction designs:

 - for our ordinal interaction (with means of mean_-1,-1 = 0, mean_-1,1 = 0, mean_1,-1 = 0, mean_1,1 = 20):  

mean_total = (0+0+0+20)/4 = 5  
mean_a=-1 = 0  
mean_a=1 = (0+20)/2 = 10  
mean_b=-1 = 0  
mean_b=1 = (0+20)/2 = 10  

sqrt(sum((mean_ab - mean_a - mean_b + 5)^2) / (2*2))  
= sqrt((0 - 0 - 0 + 5)^2 + (0 - 10 - 0 + 5)^2 + (0 - 0 - 10 + 5)^2 + (10 - 10 - 10 + 5)^2 / (2*2))  
= sqrt((4*5^2) / 4) = sqrt(5^2) = 5 = b3  


- for our hybrid interaction (with means of mean_-1,-1 = 0, mean_-1,1 = 0, mean_1,-1 = -20, mean_1,1 = 20):  

mean_total = (0 -20 + 0 + 20)/4 = 0  
mean_a=-1 = (0+0)/2 = 0  
mean_a=1 = (-20+20)/2 = 0  
mean_b=-1 = (0-20)/2 = -10  
mean_b=1 = (0+20)/2 = 10  

sqrt(sum((mean_ab - mean_a - mean_b + 0)^2) / (2*2))  
= sqrt((0 - 0 - (-10) + 0)^2 + (-20 - 0 - (-10) + 0)^2 + (0 - 0 - 10 + 0)^2 + (20 - 0 - 10 + 0)^2 / (2*2))  
= sqrt((4*10^2) / 4) = sqrt(10^2) = 10 = b3  

- for our disordinal interaction (with means of mean_-1,-1 = 20, mean_-1,1 = -20, mean_1,-1 = -20, mean_1,1 = 20):  

mean_total = (10-20-20+20)/4 = 0  
mean_a=-1 = (20-20)/2 = 0  
mean_a=1 = (-20+20)/2 = 0  
mean_b=-1 = (20-20)/2 = 0  
mean_b=1 = (-20+20)/2 = 0  

sqrt(sum((mean_ab - mean_a - mean_b + 0)^2) / (2*2))  
= sqrt((20 -0 -0 + 0)^2 + (-20 -0 -0 + 0)^2 + (-20 -0 -0 + 0)^2 + (20 -0 -0 + 0)^2 / (2*2))  
= sqrt((4*20^2) / 4) = sqrt(20^2) = 20 = b3  

(see also Cohen, J. (1977). Statistical Power Analysis for the Behavioral Sciences (Revised Ed.), Hays, 1981;
Winer, 1971; Edwards, 1972 or in short  https://webpower.psychstat.org/wiki/manual/power_of_nanova)