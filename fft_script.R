
#zad 1
plotFreqSpectr <- function(X.k, fs, xlimits = c(0, fs / 2)) {
  freq <- c(0:(length(X.k) - 1)) * fs / length(X.k)
  mod <- Mod(X.k) / (length(X.k) / 2)
  # plot_ly(x = freq, y = mod, type = 'scatter', mode = 'markers',
  #         marker = list(symbol = "line-ns-open", size = 10, line = list(width = 1)))%>%
  #   layout(title = "Frequency Spectrum", xaxis = list(title = "Frequency [Hz]", range = xlimits), yaxis = list(title = "Strength", range = c(0, max(mod))),
  #          showlegend = FALSE, hovermode = "closest")
  
  plot <-
    plot_ly(
      x = freq,
      y = mod,
      type = 'scatter',
      mode = 'markers'
    ) %>%
    layout(
      title = "Frequency Spectrum",
      xaxis = list(title = "Frequency [Hz]", range = xlimits),
      yaxis = list(title = "Strength", range = c(0, 1.5 * max(mod))),
      showlegend = FALSE,
      hovermode = "closest"
    )
  
  for (i in seq_along(freq)) {
    plot <-
      plot %>% add_trace(
        x = c(freq[i], freq[i]),
        y = c(0, mod[i]),
        type = 'scatter',
        mode = 'lines',
        line = list(color = 'red', width = 1)
      )
  }
  return(plot)
}

freq <- 4 #whole cycles per second
fs <- 10 * freq
period <- 3 / freq
tt <- seq(from = 0,
          to = period - 1 / fs,
          by = 1 / fs) #time in sec

myData <- sin(2 * pi * freq * tt)
fftMyData <- fft(myData)
plotFreqSpectr(fftMyData, fs)


freq <- 7 #whole cycles per second
fs <- 100 * freq
period <- 2 / freq
tt <- seq(from = 0,
          to = period - 1 / fs,
          by = 1 / fs) #time in sec

myData <- sawtooth(2 * pi * freq * tt, 1 / 2)
fftMyData <- fft(myData)
plotFreqSpectr(fftMyData, fs)


plot_ly(x = tt,
        y = myData,
        type = 'scatter',
        mode = 'lines') %>%
  layout(
    xaxis = list(title = ""),
    yaxis = list(title = ""),
    title = "Triangle wave"
  )


freq <- 7 #whole cycles per second
fs <- 30 * freq
period <- 2 / freq
tt <- seq(from = 0,
          to = period - 1 / fs,
          by = 1 / fs) #time in sec

myData <- square(2 * pi * freq * tt)
fftMyData <- fft(myData)
plotFreqSpectr(fftMyData, fs)

plot_ly(x = tt,
        y = myData,
        type = 'scatter',
        mode = 'lines') %>%
  layout(
    xaxis = list(title = ""),
    yaxis = list(title = ""),
    title = "Square wave"
  )


freq <- 7 #whole cycles per second
fs <- 30 * freq
period <- 2 / freq
tt <- seq(from = 0,
          to = period - 1 / fs,
          by = 1 / fs) #time in sec

myData <- rnorm(length(tt))
fftMyData <- fft(myData)
plotFreqSpectr(fftMyData, fs)

plot_ly(x = tt,
        y = myData,
        type = 'scatter',
        mode = 'lines') %>%
  layout(
    xaxis = list(title = ""),
    yaxis = list(title = ""),
    title = "Random noise"
  )


#parseval
checkPerseval <- function(sigData, fftData) {
  tolerance <- 1e-04
  sumInTime <- sum(abs(sigData) ^ 2)
  print(sumInTime)
  sumInFreq <- sum(abs(fftData) ^ 2) / length(fftData)
  print(sumInFreq)
  abs(sumInFreq - sumInTime) < tolerance
}
freq <- 3.5 #whole cycles per second
fs <- 5 * freq
period <- 2 / freq
tt <- seq(from = 0,
          to = period - 1 / fs,
          by = 1 / fs) #time in sec

signal <- sin(2 * pi * freq * tt)
signal <- signal / max(abs(signal))
signal <- signal * period
signalFft <- fft(signal)
checkPerseval(signal, signalFft)

#tricks
freq1 <- 3.05 #whole cycles per second
freq2 <- 3.5
phi1 <- pi / 2
phi2 <- pi / 4
fs <- 10 * freq1

period <- 5
timeOriginal <- seq(from = 0,
          to = period - 1 / fs,
          by = 1 / fs) #time in sec

origSignal <- sin(2 * pi * freq1 * timeOriginal + phi1) + sin(2 * pi * freq2 * timeOriginal + phi2)
plotSignal <- function(time, values) {
  plot <- plot_ly(x = time,
          y = values,
          type = 'scatter',
          mode = 'lines') %>%
    layout(
      xaxis = list(title = ""),
      yaxis = list(title = ""),
      title = ""
    )
  plot
}
# fft
origFft <- fft(origSignal)

plotSignal(timeOriginal, origSignal)

# original signal's spectra
plotFreqSpectr(origFft, fs)

# mirroring
mirroredSignal <- c(origSignal, rev(origSignal))
plotSignal(timeMirrored, mirroredSignal)

mirroredFft <- fft(mirroredSignal)
plotFreqSpectr(mirroredFft, fs)

#repeating 
repeatedSignal <-
  rep(origSignal, times = 5)
repeatedFft <- fft(repeatedSignal)
plotFreqSpectr(repeatedFft, fs)

#padding 
paddedSignal <-
  c(origSignal, rep(0, length(origSignal) * 5))
paddedFft <- fft(paddedSignal)
plotFreqSpectr(paddedFft, fs)

# library(pracma)
#real data
myData <- read.csv("C:/Users/48503/Downloads/AirPassengers.csv")
detrendedData <- detrend(myData$X.Passengers)
plot <- plot_ly(data, x = ~Month, y = ~X.Passengers, type = "scatter", mode = "lines", name = "Air Passenger Data") %>%
  layout(title = "Air Passenger Data",
         xaxis = list(title = "Month"),
         yaxis = list(title = "Number"))
plot



plot <- plot_ly(data, x = ~Month) %>%
  add_trace(y = ~X.Passengers, type = "scatter", mode = "lines", name = "Original") %>%
  add_trace(y = detrendedData, type = "scatter", mode = "lines", name = "Detrended") %>%
  layout(title = "",
         xaxis = list(title = "Month"),
         yaxis = list(title = "Number"))

plot


diffValues <- diff(detrendedData)
time <- as.numeric(as.yearmon(myData$Month))
diffTime = diff(time)
deriv <- diffValues/diffTime
plotSignal()
fftDeriv <- fft(deriv)
plotFreqSpectr(fftDeriv, 1/0.08333333)

