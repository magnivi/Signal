#include <pybind11/pybind11.h>
#include <pybind11/stl.h> 
#include <matplot/matplot.h>
#include <vector>
#include <string>
#include "AudioFile.h"

#define STRINGIFY(x) 
#define MACRO_STRINGIFY(x) STRINGIFY(x)

namespace plt = matplot;
namespace py = pybind11;
const int M_PI = 3.14;
int getNumSamples(const std::string& fileName) {
    AudioFile<short> audioFile;
    if (!audioFile.load(fileName)) {
        throw std::runtime_error("Failed to load audio file: " + fileName);
    }
    return audioFile.getNumSamplesPerChannel();
}

//ex.1 - plot audio
void process_audio(std::string fileName, int showSamples) {
    AudioFile<short> audioFile;
    audioFile.load(fileName);

    int numChannels = audioFile.getNumChannels();
    int numSamples = audioFile.getNumSamplesPerChannel();

    std::vector<short> samples;

    for (int i = 0; i < numSamples; ++i) {
        for (int channel = 0; channel < numChannels; ++channel) {
            samples.push_back(audioFile.samples[channel][i]);
        }
    }
    //if showSamples>numSamples plot all samples
    std::vector<short> subset;
    if (showSamples < numSamples) {
        for (int i = 0; i < showSamples; i++)
            subset.push_back(samples[i]);
        plt::plot(subset);
        plt::show();
    }
    else {
        plt::plot(samples);
        plt::show();
    }

}

//ex.3 - 1D Filter for vectors
//to plot all samples choose num_samples 0
void vectorFilter1D(std::vector<short>& samples, std::vector<short>& kernel1D, int num_samples) {
    std::vector<int> subset;
    if (num_samples == 0) {
        for (int i = 0; i < samples.size(); i++)
            subset.push_back(samples[i]);
    }
    else if (num_samples > 0) {
        for (int i = 0; i < num_samples; i++)
            subset.push_back(samples[i]);
    }

    std::vector<short> src;
    std::vector<short> out;
    src.push_back(0);
    for (int i : subset)
        src.push_back(i);
    src.push_back(0);

    for (int i = 0; i < subset.size(); i++) {
        int sum = 0;
        for (int j = 0; j < kernel1D.size(); j++) {
            sum += kernel1D[j] * src[i + j];
        }
        out.push_back(sum);
    }
    plt::plot(out, "*");
    plt::show();
}

//ex.3 - 1D Filter for audio
//to plot all samples choose num_samples 0
void Filter1DAudio(std::string fileName, std::vector<short>& kernel1D, int num_samples) {
    AudioFile<short> audioFile;
    audioFile.load(fileName);

    int numChannels = audioFile.getNumChannels();
    int numSamples = audioFile.getNumSamplesPerChannel();

    int readSamples = std::min(num_samples, numSamples);

    std::vector<short> samples;
    std::vector<short> src;
    std::vector<short> out;

    for (int i = 0; i < readSamples; ++i) {
        for (int channel = 0; channel < numChannels; ++channel) {
            samples.push_back(audioFile.samples[channel][i]);
        }
    }

    if (readSamples > 0) {
        std::vector<short> subset(samples.begin(), samples.begin() + readSamples);
        plt::plot(subset);
        plt::show();
    }

    src.push_back(0);
    for (int i : samples)
        src.push_back(i);
    src.push_back(0);

    for (int i = 0; i < samples.size(); i++) {
        int sum = 0;
        for (int j = 0; j < kernel1D.size(); j++) {
            sum += kernel1D[j] * src[i + j];
        }
        out.push_back(sum);
    }

    plt::plot(out);
    plt::show();
}

//ex.9 computing corelation of 2 signals
std::vector<double> calculateCorrelation(const std::vector<double>& y, const std::vector<double>& x) {
    int K = y.size() - 1;
    int M = x.size() - 1;
    int size = K + M + 1;

    std::vector<double> w(size, 0.0);

    for (int n = -K; n <= M; ++n) {
        double sum = 0.0;
        for (int k = std::max(0, -n); k <= std::min(K, M - n); ++k) {
            sum += y[k] * x[n + k];
        }
        w[n + K] = sum;
    }
    plt::plot(w);
    plt::show();
    return w;
}
void plotSignal(std::vector<double> signal) {
    plt::plot(signal);
    plt::xlabel("Samples");
    plt::ylabel("Amplitude");
    plt::title("Signal");
    plt::grid(true);
    plt::show();
}
//ex.4 - generating sin/cos/pwm.sawTooth signals
void generateWave(std::string WaveType, int frequency, double duration) {
    double samplingFrequency = 100;  //44100
    double maxDuration = 5;          // Maximum duration in seconds

    if (duration > maxDuration) {
        duration = maxDuration;
    }

    int numSamples = static_cast<int>(duration * samplingFrequency);
    std::vector<double> signal(numSamples);

    if (WaveType == "sin") {
        for (int i = 0; i < numSamples; ++i) {
            double t = static_cast<double>(i) / samplingFrequency;
            signal[i] = sin(2.0 * M_PI * frequency * t);
        }
    }
    else if (WaveType == "cos") {
        for (int i = 0; i < numSamples; ++i) {
            double t = static_cast<double>(i) / samplingFrequency;
            signal[i] = cos(2.0 * M_PI * frequency * t);
        }
    }
    else if (WaveType == "pwm") {
        double period = 1.0 / frequency;
        double halfPeriod = period / 2.0;
        for (int i = 0; i < numSamples; ++i) {
            double t = static_cast<double>(i) / samplingFrequency;
            double remainder = fmod(t, period);
            signal[i] = (remainder < halfPeriod) ? 1.0 : -1.0;
        }
    }
    else if (WaveType == "sawtooth") {
        double period = 1.0 / frequency;
        for (int i = 0; i < numSamples; ++i) {
            double t = static_cast<double>(i) / samplingFrequency;
            signal[i] = 2.0 * (t / period - floor(0.5 + t / period));
        }
    }
    else {
        std::cerr << "Unknown wave type: " << WaveType << std::endl;
        return;
    }

    plotSignal(signal);
}

PYBIND11_MODULE(Signal, m) {
    m.doc() = R"pbdoc(
        Pybind11 example plugin
        -----------------------

        .. currentmodule:: cmake_example

        .. autosummary::
           :toctree: _generate

           add
           subtract
    )pbdoc";

    m.def("vectorFilter1D", &vectorFilter1D, R"pbdoc(Plot file using matplot.)pbdoc");
    m.def("Filter1DAudio", &Filter1DAudio, R"pbdoc(Plot file using matplot.)pbdoc");
    m.def("process_audio", &process_audio, R"pbdoc(Plot file using matplot.)pbdoc");
    m.def("getNumSamples", &getNumSamples, R"pbdoc(Plot file using matplot.)pbdoc");
    m.def("calculateCorrelation", &calculateCorrelation, R"pbdoc(Plot file using matplot.)pbdoc");
    //m.def("generateWave", &generateWave, "generate wave function", py::arg("WaveType"), py::arg("frequency"), py::arg("duration"), py::arg("samplingFrequency"));
    m.def("generateWave", &generateWave, "generate wave function", py::arg("WaveType"), py::arg("frequency"), py::arg("duration"));
    //m.def("process_audio2", &process_audio, "generate wave function", py::arg("FileName"), py::arg("Samples"));



//#ifdef VERSION_INFO
//    m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
//#else
//    m.attr("__version__") = "dev";
//#endif
}
