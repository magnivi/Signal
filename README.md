

If you would like to test it in python, simply import generated .pyd file and run:

import Signal

#Testing

#print(Signal.getNumSamples("audioFile.wav"))

#Signal.process_audio("audioFile.wav", 2000)

array3 = [10,9,8,7,6,5,4,3,2,1]
kernel1D = [0,1,0]
#Signal.vectorFilter1D(array3, kernel1D, 10)

#Signal.Filter1DAudio("audioFile.wav", kernel1D, 1000)

array1 = [1,2,3,4,5,6,7,8,9,10]
array2 = [10,9,8,7,6,5,4,3,2,1]
#Signal.calculateCorrelation(array1, array2)

#waveType(sin/cos/pwm/sawTooth)  /  freqency  /  duration
#Signal.generateWave("sin", 1, 5)






