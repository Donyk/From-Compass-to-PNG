#These lines are here to import the necessary modules
print("Hello. Your favorite Simple Western online tool (and more) is now available here : https://proteinsimple-non-official.benjaminps.repl.co/ \n\n")
import os
os.environ['MPLCONFIGDIR'] = os.getcwd() + "/configs/"
import matplotlib.pyplot as plt
import matplotlib as mp
import matplotlib.gridspec as gridspec
from matplotlib import cm
from scipy.optimize import curve_fit
import numpy as np
im = plt.imread('logo.png')
from matplotlib.gridspec import GridSpec
import requests


#This defines the EC90tion which will be used to fit all the data points from the users=>
def EC90(x, a, b, c):
    global func_label
    func_label='Fit:' r'$y = b+(a-b)(\frac{1}{1+(\frac{x}{C})})$'
    global get_bounds
    def get_bounds():
      global bounds
      bounds = ((-np.inf, -np.inf, 0), (np.inf, np.inf, 1))
    global funcline
    def funcline(x):
      return (popt[2]*((popt[0]-popt[1])/((popt[1]*x)-popt[1])-1))
    return b + (a - b) / (1 + (x / c))

def PS(x, a, b, c):
    global func_label
    func_label='Fit:' r' $y = b+(a-b)(2^{-(\frac{x}{C})})$'
    global get_bounds
    def get_bounds():
      global bounds
      bounds = ((-np.inf, -np.inf, 0), (np.inf, np.inf, 1))
    global funcline
    def funcline(x):
      return (np.log(((popt[1]*x)-popt[1])/(popt[0]-popt[1]))/-np.log(2))*popt[2]
    return b + (a - b) * (2 ** -(x / c))

def KD(x, a, b, c):
    global func_label
    func_label='Fit:' r' $y = {(\frac{x}{c+x})}b$'
    global get_bounds
    def get_bounds():
      global bounds
      bounds = ((-0, 100, 0), (0.0000000001, np.inf, 1))
    global funcline
    def funcline(x):
      return popt[2]/((popt[1]/(x*popt[1]))-1)
    return a + (b * (x / (c + x)))


#Defines the lists of variables. They will be filled with inputs from the users
Columns = []
Titles = input('\nCopy/Paste all data from Compass.\nPress ENTER twice when done.\n\n')
Columns = [x.strip() for x in Titles.split('	')]
params = {"name": "benji"}
response = requests.get('https://serverless-test-nfabredev.vercel.app/api', params=params)
json = response.json()
index_name = Columns.index('Name')
index_sample = Columns.index('Sample')
index_primary = Columns.index('Primary')
index_PeakHeight = Columns.index('Height')
index_PeakArea = Columns.index('Area')
index_Baseline = Columns.index('Baseline')
index_S_N = Columns.index('S/N')
Columns = [[x] for x in Columns]
Sample = []
Primary = []
PeakHeight = []
PeakArea = []
S_N = []
Baseline = []

#This is a loop which will fill the lists with the inputs from the users=>
while True:
    candidates = input('\nData Peak ' + str(1 + len(Primary)) + ': ')
    if not candidates:
        break
    try:
        Columns = [x.strip() for x in candidates.split('	')]
        if Columns[index_name] == "":
          Fig_Title = (Columns[index_sample])
        else:
          Fig_Title = (Columns[index_sample] + " + " + Columns[index_name])
        Columns[index_primary] = Columns[index_primary].replace(',','.')
        Columns[index_PeakHeight] = Columns[index_PeakHeight].replace(',','.')
        Columns[index_PeakArea] = Columns[index_PeakArea].replace(',','.')
        Columns[index_S_N] = Columns[index_S_N].replace(',','.')
        Columns[index_Baseline] = Columns[index_Baseline].replace(',','.')
        Primary.append(Columns[index_primary])
        PeakHeight.append(float(Columns[index_PeakHeight]))
        PeakArea.append(float(Columns[index_PeakArea]))
        S_N.append(float(Columns[index_S_N]))
        Baseline.append(float(Columns[index_Baseline]))
    except ValueError:
        print('Hum, sorry, could not handle that input')

#Formats the Primary antibody to get a list of dilution factors 
if Primary[0].startswith('1:'):
  for x in range (len(Primary)):
    Primary[x] = Primary[x][2:]
for x in range (len(Primary)):
  Primary[x] = float(Primary[x])
xData = [1/x for x in Primary]


#Divides the Baseline by Peak Height for each datapoint. 
y2Data = [((x*100)/ph) for x,ph in zip(Baseline, PeakHeight)]
y2Data = [round(x,6) for x in y2Data]

#sorts the data in increasing order (if they were not already)
if xData[0] > xData[1]:
  xData.sort()
  y2Data_increasing = y2Data[::-1]
  PeakArea_increasing = PeakArea[::-1]
  S_N_increasing = S_N[::-1]
else:
   y2Data_increasing = y2Data
   PeakArea_increasing = PeakArea
   S_N_increasing = S_N

#prints (Peak Area 2 / Peak Area 1) and (Dilution Factor 2 / Dilution Factor 1)
for i in range (len(xData)-1):
  print("\nPeak Area 1:" + str(int(1/xData[i+1])) + " / Peak Area 1:" + str(int(1/xData[i])) + " = " + str(round(PeakArea_increasing[i+1]/PeakArea_increasing[i],3)))
  print("1:" + str(int(1/xData[i+1])) + " / 1:" + str(int(1/xData[i])) + " = " + str(xData[i+1]/xData[i]))

#Counts the lenth of the dataset
NbData = len(xData)
#e = np.repeat(1, NbData)


def CurvFit_Plot(func):
  global popt
  func(10, 0, 100, 0.1)
  function_name = func.__name__
  get_bounds()
  print("\n\n\n\n\n\nPreparing your Graph " + str(function_name))
  try:
    popt, pcov = curve_fit(func, xData, PeakArea_increasing, bounds=bounds)
  except RuntimeError:
    print("Fit impossible with the " + str(function_name) + "function")
    fig = plt.figure(figsize=(12,7))
    plt.text(0.1, 0.5, "Impossible fit", fontsize=72)
    fig.savefig('Graph ' + str(function_name) + '.png')
    return
  yData_increasing = [((x*100)/popt[1]) for x in PeakArea_increasing]
  yData_increasing = [round(x,3) for x in yData_increasing]
  Graph_a = round((popt[0]*100)/popt[1],3)
  Warning1font = {'family': 'serif',
          'color':  'darkred',
          'weight': 'normal',
          'size': 16,
          }
  Warning2font = {'family': 'serif',
          'color':  'darkred',
          'weight': 'normal',
          'size': 32,
          }
  Axisfont = {'color':  'k',
          'weight': 'normal',
          'size': 20,
          }
  Textfont = {'color':  'k',
          'weight': 'normal',
          'size': 12,
          }
  cbarfont = {'color':  'b',
          'weight': 'normal',
          'size': 14,
          }
  #Creates the figure
  fig = plt.figure(figsize=(12,7))
  gs1 = GridSpec(3, 7, left=0.125, right=0.975, top=0.875, bottom=0.22, wspace=0.01)
  ax1 = fig.add_subplot(gs1[:-1, :-1])
  ax2 = fig.add_subplot(gs1[-1, :])
  #The antibody titration datapoints are ploted on a graph with red dots
  ax1.plot(xData, yData_increasing, 'ro', label='Experimental data')
  #This defines the x-axis and y-axis scale for the final figure as well as the ticks and labels
  ax1.set_ylabel('Antibody Saturation\n$\\regular_{(\%\ of\ Target\ Saturated)}$', fontdict=Axisfont)
  ax1.tick_params(axis="y", labelsize=14)
  ax1.set_xticks([])
  #prints the variables a, b, c
  print("\nReal variables: ")
  print("a = " + str(int(popt[0])))
  print("b = " + str(int(popt[1])))
  print("c = " + str(round(popt[2],6)))
  #This calculates the antibody dilution needed to get 90%, 95%, 99% of saturation based on the variables (a, b, c) determined during the fit
  DotLin90 = funcline(0.90)
  DotLin95 = funcline(0.95)
  DotLin99 = funcline(0.99)
  #Plots the modeled curve on the experimental datapoints
  xFit = np.arange(0, 0.21, 0.001)
  ax1.plot(xFit, func(xFit, Graph_a, 100, popt[2]), 'r', linewidth=4, alpha=0.4, label=func_label)
  #Plots the modeled curve >90% saturation as well as the dotted vertical line
  xFit2 = np.arange(DotLin90, 0.21, 0.001)
  ax1.plot(xFit2, func(xFit2, Graph_a, 100, popt[2]), 'g', linewidth=2, label='Saturation > 90%')
  ax1.axvline(DotLin90, color='k', linewidth=2, linestyle=':')	
  #This calculates R²
  residuals = PeakArea_increasing - func(xData, *popt)
  ss_res = np.sum(residuals**2)
  ss_tot = np.sum((PeakArea_increasing - np.mean(PeakArea_increasing))**2)
  R2 = round(1 - (ss_res / ss_tot),5)
  if NbData > 3:
    ax1.plot([], [], ' ', label="$R^2$= " + str(R2))
  #Figure Title is added 
  ax1.set_title(Fig_Title, fontsize=26, loc='center', pad=20)
  if NbData <= 3:
    print("\nWarning: R² cannot be determined with 3 datapoints or fewer")
  #This calculates the antibody dilution needed to get 95% and 99% of saturation based on the variables (a, b, c) determined during the fit. Legends the graph
  if int(1/DotLin90) >= 5:
    print("\n90% antibody saturation reached at 1:" + str(int(1/DotLin90)))
  else:
    print("\nWarning: 90% saturation cannot be reached (1:" + str(round(1/DotLin90, 1)) + ")")
  if int(1/DotLin95) >= 5:
    print("95% antibody saturation reached at 1:" + str(int(1/DotLin95)))
  else:
    print("Warning: 95% saturation cannot be reached (1:" + str(round(1/DotLin95, 1)) + ")")
  if int(1/DotLin99) >= 5:
    print("99% antibody saturation reached at 1:" + str(int(1/DotLin99)))
  else:
    print("Warning: 99% saturation cannot be reached (1:" + str(round(1/DotLin99, 1)) + ")")
  #sets the limits for the x-axis for top graph
  if int(1/DotLin90) > 5 or int(1/DotLin90) < 1:
    ax1.set_xlim([0, 0.22])
  else:
    ax1.set_xlim([0, DotLin95])
  #Prints the legend on the graph 
  ax1.legend(loc=4, fontsize=16)
  #############################################################################
  #Interpolates the x and y coordinates (baseline points) in between the experimental datapoints
  xinterp = np.arange(min(xData), max(xData), 0.0002)
  yinterp = np.interp(xinterp, xData, y2Data_increasing)
  #Sets the Green-Red gradient from 5% to 20% baseline/Peak height
  norm = mp.colors.Normalize(vmin=5, vmax=20, clip=True)
  cmap = plt.cm.get_cmap('RdYlGn_r')
  my_cmap = plt.cm.RdYlGn_r(norm(yinterp))
  sm = cm.ScalarMappable(norm=norm, cmap=cmap)
  cbar = plt.colorbar(sm, ticks=[6, 19], extend='both', fraction=0.0615, pad=0.081)
  cbar.ax.set_yticklabels(['Good', 'Too \nhigh'], fontdict=cbarfont)
  #Plots the datapoints and the interpolated data
  ax2.bar(xinterp, yinterp, 0.001, color=my_cmap, alpha=1)
  ax2.plot(xData, y2Data_increasing, 'b', linewidth=1.5)
  ax2.plot(xData, y2Data_increasing, 'bx', markersize=10, label='Baseline / Peak Height (%)')
  #This defines the y-axis scale for the final figure as well as the ticks and labels
  if max(yinterp) < 20:
    ax2.set_ylim([0, 22.5])
  plt.xticks((0.005, DotLin90, 0.01, 0.02, 0.05, 0.1, 0.2), ('1:200', "≈1:" + str(int(1/DotLin90)) + "------", '1:100', '1:50', '1:20', '1:10', '1:5'), color='k', size=14, rotation=90)
  ax2.tick_params(axis="x", labelsize=14, rotation=90)
  ax2.tick_params(axis="y", labelsize=14)
  ax2.set_ylabel('Baseline\n$\\regular_{(\%\ of\ Peak\ Height)}$', fontdict=Axisfont)
  ax2.yaxis.label.set_color('blue')
  ax2.axvline(DotLin90, color='k', linewidth=2, linestyle=':')	
  plt.gca().get_xticklabels()[1].set_color('green') 
  ax2.set_xlabel('Dilution of Primary Antibody', fontdict=Axisfont)
  ax2.tick_params(axis='y', colors='blue')
  #Plots the Signal/Noise
  ax3 = ax2.twinx()
  ax3.set_ylim([0, 1.5*max(S_N_increasing)])
  ax3.plot(xData, S_N_increasing, 'rx', markersize=10, label='Signal / Noise')
  ax3.set_ylabel('Signal / Noise', fontdict=Textfont)
  ax3.yaxis.label.set_color('red')
  ax3.tick_params(axis='y', colors='red', labelsize=14)
  ax3.tick_params(axis="y",direction="in", pad=-40)
  ax3.yaxis.label.set_color('red')
  #organise the Second y-axis
  if max(S_N_increasing) > 10000:
    Ben = round(max(S_N_increasing), -4)
  elif max(S_N_increasing) > 1000:
    Ben = round(max(S_N_increasing), -3)
  elif max(S_N_increasing) >100:
    Ben = round(max(S_N_increasing), -2)
  else:
    Ben = round(max(S_N_increasing), -1)
  ax3.set_yticks([Ben, Ben*0.75, Ben*0.5, Ben*0.25])
  ax2.legend(loc=2, fontsize=12)
  ax3.legend(loc=1, fontsize=12)
  #sets the limits for the x-axis for bottom graph
  if int(1/DotLin90) > 5 or int(1/DotLin90) < 1:
    ax2.set_xlim([0, 0.22])
  else:
    ax2.set_xlim([0, DotLin95])
  xlim = ax1.get_xlim()
  xlim1percent = (xlim[1]+(xlim[1]*0.01))
  #Change the Y-axis scale if the saturation cannot be reached
  if int(1/DotLin90) < 1: # and R2 < 0.9:
    ylim = ax1.get_ylim()
    #yticks = np.arange(0, max(yData_increasing), 10)
    #ax1.set_yticks(yticks)
    ax1.tick_params(axis='y', colors='darkred', labelsize=18)
    ylim_scale = (ylim[1]-ylim[0])
    if R2 < 0.9 or (1/DotLin90) < 0.5:
      ax1.text((xlim[0]+(xlim[1]*0.01)), (ylim[0]+(0.85*ylim_scale)),"Impossible fit", fontdict=Warning2font)
    else:
      ax1.text((xlim[0]+(xlim[1]*0.01)), (ylim[0]+(0.85*ylim_scale)),"Saturation impossible", fontdict=Warning2font)
      print("\nSaturation impossible \nPlease consider the following: \n- Diluting the sample \n- Increasing the antibody incubation time (120 minutes maximum) \n- Using an alternative antibody \n- Contacting your F.A.S. if you need support\n")
  else:
    ax1.set_ylim([0, 102.5])
    ylim = ax1.get_ylim()
    ylim_scale = (ylim[1]-ylim[0])
  #Labels the graph (on the right side)
  if int(max(PeakHeight)) > 300000:
    print("\nWARNING: Risks of burnouts detected. Please check all-exposures or contact your FAS.")
    ax1.text((xlim[0]+(xlim[1]*0.01)), (ylim[0]-(0.075*ylim_scale)),"WARNING: Risks of burnouts detected. Please check all-exposures.", fontdict=Warning1font)
  ax1.text((1.01*(xlim[1]-xlim[0])), ylim[0]+0.85*ylim_scale,"Graph variables: ", fontdict=Textfont)
  ax1.text(xlim1percent, ylim[0]+0.75*ylim_scale,"a = " + str(round(Graph_a,4)), fontdict=Textfont)
  ax1.text(xlim1percent, ylim[0]+0.65*ylim_scale,"b = 100", fontdict=Textfont)
  ax1.text(xlim1percent, ylim[0]+0.55*ylim_scale,"c = " + str(round(popt[2],4)), fontdict=Textfont)
  if int(1/DotLin90) >= 5:
    ax1.text(xlim1percent, ylim[0]+0.35*ylim_scale,"1:" + str(int(1/DotLin90)) + " → 90%", fontdict=Textfont)
  else:
    ax1.text(xlim1percent, ylim[0]+0.35*ylim_scale,"90% not reachable", fontdict=Textfont)
  if int(1/DotLin95) >= 5:
    ax1.text(xlim1percent, ylim[0]+0.25*ylim_scale,"1:" + str(int(1/DotLin95)) + " → 95%", fontdict=Textfont)
  else:
    ax1.text(xlim1percent, ylim[0]+0.25*ylim_scale,"95% not reachable", fontdict=Textfont)
  if int(1/DotLin99) >= 5:
    ax1.text(xlim1percent, ylim[0]+0.15*ylim_scale,"1:" + str(int(1/DotLin99)) + " → 99%", fontdict=Textfont)
  else:
    ax1.text(xlim1percent, ylim[0]+0.15*ylim_scale,"99% not reachable", fontdict=Textfont)
  ax1.text((xlim[1]+(xlim[1]*0.07)), ylim[0]-0.525*ylim_scale,"Baseline scale", fontdict=cbarfont, rotation=90)
  #It creates the final figure and saves it
  fig.figimage(im, 300, 320, alpha=.2)
  fig.savefig('Graph ' + str(function_name) + '.png')
  print("\nYour Graph " + str(function_name) + " is ready.")

#CurvFit_Plot(EC90)
CurvFit_Plot(PS)
#CurvFit_Plot(KD)


print("\n\n\n\n\nIf you like this tool, feedback is greatly appreaciated : benjamin.djian@bio-techne.com\nHave a great day!")

plt.show()



