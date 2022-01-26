import numpy as np
import math
import matplotlib.pyplot as plt

# given a file_name without file extension
# extract the detector dimensions form file_name.att
# import and reshape the detector singal from file_name.bin
# return the dimensions of the detector and the detector signal
def import_detector_signal_from_file(file_name):
    attribute_file_name = file_name + '.att'
    attribute_file = open(attribute_file_name, "r")
    line = attribute_file.readline().strip()
    while line:
        split_line = line.partition("=")
        if(split_line[0] == 'Pixels total signal [Nx, Ny] '):
            lines_removed_brackets = split_line[2].strip('[] ')
            values = lines_removed_brackets.partition(',')
            width = int(values[0])
            height = int(values[2])
        line=attribute_file.readline().strip()
    attribute_file.close()
    binary_file_name = file_name + '.bin'
    Detector_Signal = np.reshape(np.fromfile(binary_file_name,float),(height,width))
    return [width, height, Detector_Signal]



def get_detector_dimensions(file_name):
    attribute_file_name = file_name + '.att'
    attribute_file = open(attribute_file_name, "r")
    line = attribute_file.readline().strip()
    print(line)
    while line:
        print(line)
        split_line = line.partition("=")
        if(split_line[0] == 'Dimensions [x,y] '):
            lines_removed_brackets = split_line[2].strip('[] ')
            values = lines_removed_brackets.partition(',')
            FOV_x = float(values[0])
            FOV_y = float(values[2])
            break
        line = attribute_file.readline().strip()
    attribute_file.close()
    return [FOV_x, FOV_y]

def get_number_of_histories(file_name):
    attribute_file_name = file_name + '.att'
    attribute_file = open(attribute_file_name, "r")
    line = attribute_file.readline().strip()
    while line:
        split_line = line.partition("=")
        if(split_line[0] == 'Number of histories '):
            number_of_histories = int(split_line[2])
            break
        line = attribute_file.readline().strip()
    attribute_file.close()
    return number_of_histories

#function computes phase stepping curves for a binary analyzer grating with duty cycle 0.5
# using the following parameters:
#   - nPixelX, nPixelY; number of pixels in x and y directions for rebinning the input signal
#       IMPORTANT: array dimensions are assumed to be an integer multiple of nPixelX and nPixelY
#   - nMask; the width of the mask, i.e. whow many consecutive MC data points are covered by the grating
#       IMPORTANT: the width of the grating sections is assumed to be an integer multiple of the width of a MC data point
#   - nPhaseSteps number of phase steps to perform
#   - Detector_Signal input signal, which can be generated by import_detector_signal_from_file
#! IMPORTANT: Does NOT work with 1 dimensional Detector_Signal input.
def calculate_phase_stepping_curve(nPixelX, nPixelY, nMask, nPhaseSteps, Detector_Signal):
    PhaseSteppingCurve = np.zeros([nPhaseSteps, nPixelY, nPixelX])

    #get dimensions of the detector signal
    width = Detector_Signal.shape[1]
    height = Detector_Signal.shape[0]

    #How many MC data points make a pixel
    nDataPointsPerPixelX = math.floor(width/nPixelX)
    nDataPointsPerPixelY = math.floor(height/nPixelY)

    #Period of G2 in units of MC pixels
    nPeriodG2 = 2 * nMask
    #how many MC-pixels to shift the grating each phase step
    nOffset=math.floor(2*nMask/nPhaseSteps)
    print("nOffset: " + str(nOffset))
    #Calculate the phase stepping curve for each pixel called MacroPixel
    for nPixelCounterX in range(nPixelX):
        for nPixelCounterY in range(nPixelY):
            MacroPixel = np.copy(Detector_Signal[(nPixelCounterY * nDataPointsPerPixelY) : ((nPixelCounterY+1) * nDataPointsPerPixelY) , (nPixelCounterX * nDataPointsPerPixelX) : ((nPixelCounterX+1) * nDataPointsPerPixelX)])
            for nPhaseSetpCounter in range(nPhaseSteps):
                MaskedPixel = np.copy(MacroPixel)
                for nPeriodCounter in range(math.ceil(nDataPointsPerPixelX / nPeriodG2) + 2):
                    #mask pixel from 'nfrom' to 'nto'
                    nfrom = nPeriodCounter * nPeriodG2 - (nDataPointsPerPixelX * nPixelCounterX - nOffset * nPhaseSetpCounter) % nPeriodG2
                    nto = nfrom + nMask
                    if nto > nDataPointsPerPixelX:
                        nto = nDataPointsPerPixelX
                    elif nto <= 0:
                        continue

                    if nfrom <= 0:
                        nfrom = 0

                    if nfrom < nDataPointsPerPixelX:
                        MaskedPixel[:,nfrom:nto] = 0;

                PhaseSteppingCurve[nPhaseSetpCounter, nPixelCounterY, nPixelCounterX] = np.sum(np.sum(MaskedPixel))
    del MacroPixel
    del MaskedPixel
    return PhaseSteppingCurve;



#given a phase stepping curve (e.g. generated by calculate_phase_stepping_curve), get the the visibility map
#v = (Imax-Imin)/(Imax+Imin)
#def get_visibility_map(phase_stepping_curve):
#    nPixelX = phase_stepping_curve.shape[2]
#    nPixelY = phase_stepping_curve.shape[1]
#    visibility_map = np.zeros([nPixelY, nPixelX])
#    for nPixelCounterX in range(nPixelX):
#        for nPixelCounterY in range(nPixelY):
#            MaxSiganl = np.max(phase_stepping_curve[:,nPixelCounterY, nPixelCounterX])
#            MinSignal = np.min(phase_stepping_curve[:,nPixelCounterY, nPixelCounterX])
#            visibility_map[nPixelCounterY,nPixelCounterX] = (MaxSiganl - MinSignal) / (MaxSiganl + MinSignal)
#    return visibility_map


def get_visibility_map(phase_stepping_curve):
    MaxSiganl = np.max(phase_stepping_curve,axis=0)
    MinSignal = np.min(phase_stepping_curve, axis=0)
    return (MaxSignal - MinSignal) /(MaxSiganl + MinSignal)






#a function to extract the simulation time from the MC output.
#'MC_output_file_name' has to include the file extension.
def get_simulation_time(MC_output_file_name):
    cpu_time = 0.0;
    elapsed_time = 0.0;
    MC_output_file = open(MC_output_file_name, "r")
    for line in MC_output_file:
        line = line.strip();
        split_line = line.partition(":")
        if(split_line[0] == 'CPU time'):
            a = split_line[2].partition('s')
            value = a[0].strip(' ')
            cpu_time = float(value)
        elif(split_line[0] == 'Elapsed time'):
            a = split_line[2].partition('s')
            value = a[0].strip(' ')
            elapsed_time = float(value)
    MC_output_file.close()
    return [cpu_time, elapsed_time]



#a fct for determining the period of the signal
#as input assume a one dimensional object
def get_periodicity__(Averaged_Detector_Signal, d_in):
    fft_detector_signal = np.fft.fft(Averaged_Detector_Signal)
    #print(fft_detector_signal)
    n = Averaged_Detector_Signal.size
    frequencies = np.fft.fftfreq(n,d=d_in)
    #print(frequencies)
    plt.plot(frequencies, fft_detector_signal.real)
    plt.plot(frequencies, fft_detector_signal.imag)
    plt.show()
    max_index = np.argmax(np.abs(fft_detector_signal[1:-1]))
    print("max_signal: " + str(max_index))
    print("f[may_index]: " + str(fft_detector_signal[max_index + 1]))
    print("frequency: " + str(frequencies[max_index + 1]))
    print("period: " + str(1.0 / frequencies[max_index + 1]))
    return abs(1.0 / frequencies[max_index + 1])


def get_periodicity2__(Averaged_Detector_Signal, d_in):
    fft_detector_signal = np.fft.fft(Averaged_Detector_Signal)
    #print(fft_detector_signal)
    n = Averaged_Detector_Signal.size
    frequencies = np.fft.fftfreq(n,d=d_in)
    erase = int(round(n/100))
    fft_detector_signal[0:erase] = 0.0
    #print(frequencies)
    #plt.plot(frequencies, fft_detector_signal.real)
    #plt.plot(frequencies, fft_detector_signal.imag)
    #plt.show()
    max_index = np.argmax(np.abs(fft_detector_signal[1:-1]))
    freq = frequencies[max_index + 1]
    print("max_signal: " + str(max_index))
    print("f[may_index]: " + str(fft_detector_signal[max_index + 1]))
    print("frequency: " + str(frequencies[max_index + 1]))
    print("period: " + str(1.0 / frequencies[max_index + 1]))
    return abs(1.0 / frequencies[max_index + 1])
