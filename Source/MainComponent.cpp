#include "AudioFile.h"
#include "MainComponent.h"
#include "filt.h"

#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <map>
#include <stdio.h>
#include <iomanip>
#include <cmath>
#define _USE_MATH_DEFINES
#include <math.h>
#include <random>
#include <regex>
#include <filesystem>
#define PI 3.14159265
using namespace std;




class Ray { 
public:
    float x; float y; float z; float rayDirectionX; float rayDirectionY; float rayDirectionZ; double numberOfBounces = 0; double iterationsTravelled;
    vector<double> sample;
};
class ThreeDimensionalRoom {
public:
    int width; int height; int length; double wallDecayFactor;
};


class microphone {
public:
    int micLocX; int micLocY; int micLocZ; int size;
    AudioFile<double>::AudioBuffer buffer; vector<Ray> RaysDetected; int whereInBuffer;
    int micType; float micAzimuthalAngle; float micPolarAngle; float micPolarAngle180Degrees; float micAzimuthalAngle180Degrees;
    float azimuthalAngle;  float polarAngle;

    void raysDetectedByMic(Ray ray) {
        micType = 0;

        /*Default to Omni*/
        switch (micType) {
        case 0:
            omni(ray);
            break;
        case 1:
            cardioid(ray);
            break;
        case 2:
            figure8(ray);
            break;
        }
    }
    void omni(Ray ray) {
        RaysDetected.push_back(ray);
    }

    void cardioid(Ray ray) {
        azimuthalAngle = 0; polarAngle = 0;

        azimuthalAngle = getAzimuthal(ray); polarAngle = getPolar(ray);


        //DBG("Incoming Ray Azumuthal Angle");
        //DBG(azimuthalAngle);
        //DBG("Incoming Ray Polar Angle");
        //DBG(polarAngle);

        micPolarAngle = 15; micAzimuthalAngle = 15;


        /* This solves ^ problem */
        if ((polarAngle > (micPolarAngle - 75)) && (polarAngle < micPolarAngle + 75)) {
            if ((azimuthalAngle > (micAzimuthalAngle - 75)) && (azimuthalAngle < (micAzimuthalAngle + 75))) {
                //DBG("Mic active region therefore no transformation is needed");
                RaysDetected.push_back(ray);
            }
        }
        else if ((polarAngle < (micPolarAngle - 75)) && (polarAngle > micPolarAngle + 75) && (polarAngle > (micPolarAngle - 120)) && (polarAngle < micPolarAngle + 120)) {
            if ((azimuthalAngle > (micAzimuthalAngle - 75)) && (azimuthalAngle < (micAzimuthalAngle + 75)) && (azimuthalAngle > (micAzimuthalAngle - 120)) && (azimuthalAngle < (micAzimuthalAngle + 120))) {
                //DBG("Mic active region dimished");
                for (int j = 0; j < ray.sample.size(); j++) {
                    ray.sample.at(j) = ray.sample.at(j) / 3;
                }
                RaysDetected.push_back(ray);
            }
        }
        else {
            for (int j = 0; j < ray.sample.size(); j++) {
                ray.sample.at(j) = -ray.sample.at(j);
            }
            RaysDetected.push_back(ray);
        }
    }
    void figure8(Ray ray) {
        azimuthalAngle = 0; polarAngle = 0;
        azimuthalAngle = getAzimuthal(ray); polarAngle = getPolar(ray);


        //DBG("Incoming Ray Azumuthal Angle");
        //DBG(azimuthalAngle);
        //DBG("Incoming Ray Polar Angle");
        //DBG(polarAngle);

        micPolarAngle = 15; micAzimuthalAngle = 15;

        float micPolarAngleTemp; float micAzimuthalAngleTemp;

        micPolarAngleTemp = micPolarAngle - 180;
        micAzimuthalAngleTemp - micAzimuthalAngle - 180;
        micPolarAngleTemp = normalize(micPolarAngleTemp, 0, 360);
        micAzimuthalAngleTemp = normalize(micAzimuthalAngleTemp, 0, 360);

        if ((polarAngle > (micPolarAngle - 75)) && (polarAngle < micPolarAngle + 75)) {
            if ((azimuthalAngle > (micAzimuthalAngle - 75)) && (azimuthalAngle < (micAzimuthalAngle + 75))) {
                //DBG("Mic active Side 1 of Figure 8");
                RaysDetected.push_back(ray);
            }
        } /* This should create the other side of the figure 8 microphone. */
        /* Here we need to create a transform  for the alternative angle */
        else if ((polarAngle > (micPolarAngleTemp - 75)) && (polarAngle < micPolarAngleTemp + 75)) {
            if ((azimuthalAngle > (micAzimuthalAngleTemp - 75)) && (azimuthalAngle < (micAzimuthalAngleTemp + 75))) {
                //DBG("Mic active Side 2 of Figure 8");
                RaysDetected.push_back(ray);
            }
        } else {
            for (int j = 0; j < ray.sample.size(); j++) {
                ray.sample.at(j) = ray.sample.at(j) / 5;
            }
            RaysDetected.push_back(ray);
        }
    }

    double applyAmplitudeCalculationToSample(double sample, int numberOfBounces, int iterationsTravelled, int reverbAmount) {
        /* http://hyperphysics.phy-astr.gsu.edu/hbase/Sound/reflec.html */
        if ((double(numberOfBounces % 2)) == 1) {
            sample = -sample;
        }
        if (numberOfBounces >= 1) { 
           sample = sample / pow(2, numberOfBounces);

        }
        return sample;
    }
    /* this should return a number of rays in a vector of doubles.*/
    /* This method needs to process microphone input through rays detected and write the output audio buffer */
    void sumOfEveryBitOfDataHittingMicAtThisIteration(int currentIteration, int numSamplesPerRay, int numSamplesInChannel, int reverbAmount) {
        /* you need to sum all samples for every iteration. This is vector of vectors of doubles */
        int numRaysDetected = RaysDetected.size();
        vector<double> distanceEachRayHasTravelled;
        vector<double> numOfBouncesEachRay;
        vector<double> OutputSamplesToWriteToBuffer;

        for (int k = 0; k < numSamplesPerRay; k++) {
            OutputSamplesToWriteToBuffer.push_back(0);
        }
        for (int i = 0; i < numRaysDetected; i++) {
            distanceEachRayHasTravelled.push_back(RaysDetected[i].iterationsTravelled);
            numOfBouncesEachRay.push_back(RaysDetected[i].numberOfBounces);
        }
        
        /* Create multiple filters with less severe low pass for more realistic emulation of audio such that it low pass increases with number of bounces */
        Filter* my_filter;
        my_filter = new Filter(LPF, 51, 44.1, 2.0);
        for (int i = 0; i < numRaysDetected; i++) {
            for (int j = 0; j < numSamplesPerRay; j++) {
                if (RaysDetected[i].numberOfBounces == 1) {
                    /* https://www.cardinalpeak.com/blog/a-c-class-to-implement-low-pass-high-pass-and-band-pass-filters */
                    RaysDetected[i].sample[j] = my_filter->do_sample(RaysDetected[i].sample[j]);
                }
                RaysDetected[i].sample[j] = applyAmplitudeCalculationToSample(RaysDetected[i].sample[j], numOfBouncesEachRay[i], distanceEachRayHasTravelled[i], reverbAmount);
            }
        }
        for (int o = 0; o < numRaysDetected; o++) {
            for (int p = 0; p < numSamplesPerRay; p++) {
                OutputSamplesToWriteToBuffer[p] = OutputSamplesToWriteToBuffer[p] + RaysDetected[o].sample[p];
            }
        }
        for (int l = 0; l < OutputSamplesToWriteToBuffer.size(); l++) {
            whereInBuffer++;
            if (whereInBuffer < numSamplesInChannel) {
                 buffer[0][whereInBuffer] = buffer[0][whereInBuffer] + OutputSamplesToWriteToBuffer[l];
                 buffer[1][whereInBuffer] = buffer[1][whereInBuffer] + OutputSamplesToWriteToBuffer[l];
            }
            else {
                break;
            }
        }
        /* clear all the output samples */
        OutputSamplesToWriteToBuffer.clear();
        /* clear all the rays detected */
        RaysDetected.clear();
        /* clear the number of rays detected counter */
        numRaysDetected = 0;
    }
    float getAzimuthal(Ray ray) {
        float Azimuthal;
        Azimuthal = atan2(ray.rayDirectionX, ray.rayDirectionY) * 180 / PI;
        Azimuthal = Azimuthal + 180;
        return Azimuthal;
    }
    float getPolar(Ray ray) {
        float Polar;
        Polar = atan2(ray.rayDirectionX, ray.rayDirectionZ) * 180 / PI;
        Polar = Polar + 180;
        return Polar;
    }
    // function to calculate the distance between centre and the point
    // (cx, cy, cz) = center coords.
    // (x, y, z) = point coords
    float check(float cx, float cy, float cz, float x, float y, float z)
    {
        float x1 = pow((x - cx), 2);
        float y1 = pow((y - cy), 2);
        float z1 = pow((z - cz), 2);

        // distance between the centre
        // and given point
        return (x1 + y1 + z1);
    }

    /* https://stackoverflow.com/questions/1628386/normalise-orientation-between-0-and-360 */
    // Normalizes any number to an arbitrary range 
    // by assuming the range wraps around when going below min or above max 
    float normalize(const float value, const float start, const float end) {
        const float width = end - start;   // 
        const float offsetValue = value - start;   // value relative to 0

        return (offsetValue - (floor(offsetValue / width) * width)) + start;
        // + start to reset back to start of original range
    }
};
/* This is the sound source */
class soundSource {
public:
    int x; int y; int z; int iterator; AudioFile<double> input; vector<double> allInputSamples; int reverbAmount;
    vector<Ray> allRay;
    vector<microphone> allMic;
    vector<vector<double>> RaySamples;
    vector<float> directionX, directionY, directionZ;

    /* we went with the fibonacci spiral sphere https://bduvenhage.me/geometry/2019/07/31/generating-equidistant-vectors.html */
    void generateDirectionsFibonacciSpiralSphere(const int numDirectionsPerRay, int distanceRayTravelled) {
        vector<float> xFloat; float xTemp;
        vector<float> yFloat; float yTemp;
        vector<float> zFloat; float zTemp;
        const double gr = (sqrt(5.0) + 1.0) / 2.0;  // golden ratio = 1.6180339887498948482
        const double ga = (2.0 - gr) * (2.0 * M_PI);  // golden angle = 2.39996322972865332
        for (size_t i = 1; i <= numDirectionsPerRay; ++i) {
            const double lat = asin(-1.0 + 2.0 * double(i) / (numDirectionsPerRay + 1));
            const double lon = ga * i;
            const double x = cos(lon) * cos(lat);
            const double y = sin(lon) * cos(lat);
            const double z = sin(lat);
            xTemp = x * distanceRayTravelled; yTemp = y * distanceRayTravelled; zTemp = z * distanceRayTravelled;
            xFloat.push_back(xTemp); yFloat.push_back(yTemp); zFloat.push_back(zTemp);
            directionX = xFloat; directionY = yFloat; directionZ = zFloat;
        }
    }

    void performRayTrace(ThreeDimensionalRoom room, int numSamplesPerRay, int maxRaysInExistence, int numDirectionsPerRay, float distanceRayTravelled) {
        int numSamplesInChannel = input.getNumSamplesPerChannel();
        //int size = RaySamples.size();
        //const int numDirectionsPerRay = 75;
        
        //int distanceRayTravelled = 10;
        generateDirectionsFibonacciSpiralSphere(numDirectionsPerRay, distanceRayTravelled);

        /* Creates a vector of samples to form input buffer (which is just a vector of doubles) */
        for (int i = 0; i < numSamplesInChannel; i++) {
            allInputSamples.push_back(input.samples[0][i]);
            allInputSamples[i] = (allInputSamples[i] + input.samples[1][i]) / 2;
            /* This is to minimize chances of audio clipping given maximum value sample audio input */

            /* Scale factor from SKEW */
            /* More directions per ray, the more this number decreases */
            allInputSamples[i] = allInputSamples[i] * 0.05;
        }
        /* Now we need to seperate the input samples into groups of size n (You must deal with the remainder) */
        int numberOfRaysToCast = allInputSamples.size() / numSamplesPerRay;
        int remainder = allInputSamples.size() % numSamplesPerRay;

        /* This ensures all groups of size n are of size n to avoid vector subscript out of range error */
        for (int j = 0; j < remainder; j++) {
            allInputSamples.push_back(0);
        }
        /* Now we need to create the RaysToCast containing subdivisions of samples */
        for (int k = 0; k < numberOfRaysToCast; k++) {
            vector<double> oneRay;
            int l = k * numSamplesPerRay;
            for (int m = l; m < l + numSamplesPerRay; m++) {
                oneRay.push_back(allInputSamples[m]);
            }
            RaySamples.push_back(oneRay);
            oneRay.clear();
        }
        /* Now we need to cast all the rays */
        for (int n = 0; n < numberOfRaysToCast; n++) {
            castRay(room, RaySamples[n], n, numSamplesPerRay, numSamplesInChannel, numDirectionsPerRay, maxRaysInExistence);
            cout << "Current Iteration: " << n << " out of: " << numberOfRaysToCast << endl;
        }

    }
    void castRay(ThreeDimensionalRoom room, vector<double> vectorOfSamples, int currentIteration, int numSamplesPerRay, int numSamplesInChannel, int numDirectionsPerRay, int maxRaysInExistence) {
        float sizeOfRay = vectorOfSamples.size();
        for (int i = 0; i < directionX.size(); i++) {
            Ray ray;
            ray.sample = vectorOfSamples;
            ray.x = x; ray.rayDirectionX = directionX[i];
            ray.y = y; ray.rayDirectionY = directionY[i];
            ray.z = z; ray.rayDirectionZ = directionZ[i];
            allRay.push_back(ray); //cout << "allray size  " << allRay.size() << endl;
        }
        if (allRay.size() > maxRaysInExistence) {
            allRay.erase(allRay.begin(), allRay.begin() + numDirectionsPerRay);
        }
        for (signed int i = 0; i < allRay.size(); i++) {
            if (allRay.at(i).numberOfBounces == 10) {
                allRay.erase(allRay.begin() + i);
            }
            allRay.at(i) = traceOneRay(allRay.at(i), room);
        }
        for (int i = 0; i < allMic.size(); i++) {
            allMic[i].sumOfEveryBitOfDataHittingMicAtThisIteration(currentIteration, numSamplesPerRay, numSamplesInChannel, reverbAmount);
        }
        iterator++;
    }
    Ray traceOneRay(Ray ray, ThreeDimensionalRoom room) {
        if (ray.x >= room.width) {
            ray.rayDirectionX = -ray.rayDirectionX; ray.numberOfBounces += 1;
        }
        if (ray.x <= 0) {
            ray.rayDirectionX = -ray.rayDirectionX; ray.numberOfBounces += 1;
        }
        if (ray.y >= room.height) {
            ray.rayDirectionY = -ray.rayDirectionY; ray.numberOfBounces += 1;
        }
        if (ray.y <= 0) {
            ray.rayDirectionY = -ray.rayDirectionY; ray.numberOfBounces += 1;
        }
        if (ray.z >= room.length) {
            ray.rayDirectionZ = -ray.rayDirectionZ; ray.numberOfBounces += 1;
        }
        if (ray.z <= 0) {
            ray.rayDirectionZ = -ray.rayDirectionZ; ray.numberOfBounces += 1;
        }
        ray.iterationsTravelled++; ray.x = ray.x + ray.rayDirectionX; ray.y = ray.y + ray.rayDirectionY; ray.z = ray.z + ray.rayDirectionZ;
        for (auto& i : allMic) {

            float ans = 0;
            ans = i.check(i.micLocX, i.micLocY, i.micLocZ, ray.x, ray.y, ray.z);

            i.micLocX;
            i.micLocY;
            i.micLocZ;

            if (ans < (i.size * i.size)) {
                i.raysDetectedByMic(ray);
            }
        } return ray;
    }

    AudioFile<double>::AudioBuffer createStereoAudio(AudioFile<double>::AudioBuffer bufferLeftMicrophone, AudioFile<double>::AudioBuffer bufferRightMicrophone, AudioFile<double>::AudioBuffer outputBuffer) {       
        outputBuffer[0] = bufferLeftMicrophone[0];
        outputBuffer[1] = bufferRightMicrophone[1];
        return outputBuffer;
    }
};

//using namespace juce;
#undef DONT_SET_USING_JUCE_NAMESPACE
//==============================================================================
MainComponent::MainComponent() : juce::AudioAppComponent(otherDeviceManager), state(Stopped), /*stateOutput(StoppedOutput),*/ openButton("Open File"), playButton("Play"), stopButton("Stop"),
    startRayTraceButton("Start Ray Trace and Save To File.")
{

    otherDeviceManager.initialise(2, 2, nullptr, true);
    audioSettings.reset(new AudioDeviceSelectorComponent(otherDeviceManager, 0, 0, 0, 2, false , false, true, false));
    addAndMakeVisible(audioSettings.get());
    playbackSettingsLabel.setText("Audio Playback Settings", dontSendNotification);
    addAndMakeVisible(playbackSettingsLabel);

    // Make sure you set the size of the component after
    // you add any child components.
    setSize (800, 900);

    //Open Button
    openButton.onClick = [this] { openButtonClicked(); };
    addAndMakeVisible(&openButton);

    playButton.onClick = [this] { playButtonClicked(); };
    playButton.setColour(juce::TextButton::buttonColourId, juce::Colours::green);
    playButton.setEnabled(false);
    addAndMakeVisible(&playButton);

    stopButton.onClick = [this] { stopButtonClicked(); };
    stopButton.setColour(juce::TextButton::buttonColourId, juce::Colours::red);
    stopButton.setEnabled(false);
    addAndMakeVisible(&stopButton);



    getLookAndFeel().setColour(Slider::thumbColourId, Colours::blue);
    getLookAndFeel().setColour(Slider::rotarySliderOutlineColourId, Colours::black);
    getLookAndFeel().setColour(Slider::rotarySliderFillColourId, Colours::azure);


    maxGainDialSize = 500;
    gainLabel.setText("Gain", dontSendNotification);
    addAndMakeVisible(gainLabel);
    gainDial.setSliderStyle(Slider::SliderStyle::Rotary);
    gainDial.setTextBoxStyle(Slider::TextBoxBelow, true, 200, 25);
    gainDial.setRange(0, maxGainDialSize);
    gainDial.setValue(1);
    gainDial.setSkewFactorFromMidPoint(maxGainDialSize/16);
    gainDial.addListener(this);
    addAndMakeVisible(gainDial);


    /* Reverberation Chamber Size XYZ */
    maxRoomSize = 1000;
    /* Size X */
    dial1.setSliderStyle(Slider::SliderStyle::Rotary);
    dial1.setTextBoxStyle(Slider::TextBoxAbove, true, 75, 25);
    dial1.setRange(50, maxRoomSize);
    dial1.setValue(70);
    dial1.setTextValueSuffix("m");
    dial1.setSkewFactorFromMidPoint(maxRoomSize / 4);
    dial1.addListener(this);
    roomSizeXYZ.setText("Reverberation Chamber Size (Metres) X/Y/Z", dontSendNotification);

    /* Size Y */
    dial2.setSliderStyle(Slider::SliderStyle::Rotary);
    dial2.setTextBoxStyle(Slider::TextBoxAbove, true, 75, 25);
    dial2.setRange(30, maxRoomSize);
    dial2.setValue(70);
    dial2.setTextValueSuffix("m");
    dial2.setSkewFactorFromMidPoint(maxRoomSize / 4);
    dial2.addListener(this);

    /* Size Z */
    dial3.setSliderStyle(Slider::SliderStyle::Rotary);
    dial3.setTextBoxStyle(Slider::TextBoxAbove, true, 75, 25);
    dial3.setRange(50, maxRoomSize);
    dial3.setValue(70);
    dial3.setTextValueSuffix("m");
    dial3.setSkewFactorFromMidPoint(maxRoomSize / 4);
    dial3.addListener(this);

    addAndMakeVisible(roomSizeXYZ);
    addAndMakeVisible(dial1);
    addAndMakeVisible(dial2);
    addAndMakeVisible(dial3);

    soundSource sound;
    sound.x = d4_SoundSourceLocX;
    sound.y = d5_SoundSourceLocY;
    sound.z = d6_SoundSourceLocZ;
    sound.reverbAmount = 25;


    std::vector<microphone> allMics;
    microphone StereoLeft;
    microphone StereoRight;


    /* Sound Source Location */
    /* location X */
    dial4.setSliderStyle(Slider::SliderStyle::Rotary);
    dial4.setTextBoxStyle(Slider::TextBoxAbove, true, 200, 25);
    dial4.setRange(0, dial1.getValue());
    dial4.setValue(10);
    dial4.setTextValueSuffix("m");
    soundSourceLocationXYZ.setText("Sound Source Location Coordinates X/Y/Z ", dontSendNotification);
    addAndMakeVisible(soundSourceLocationXYZ);

    /* location Y */
    dial5.setSliderStyle(Slider::SliderStyle::Rotary);
    dial5.setTextBoxStyle(Slider::TextBoxAbove, true, 75, 25);
    dial5.setRange(0, dial2.getValue());
    dial5.setValue(10);
    dial5.setTextValueSuffix("m");

    /* location Z */
    dial6.setSliderStyle(Slider::SliderStyle::Rotary);
    dial6.setTextBoxStyle(Slider::TextBoxAbove, true, 75, 25);
    dial6.setRange(0, dial3.getValue());
    dial6.setValue(10);
    dial6.setTextValueSuffix("m");

    addAndMakeVisible(dial4);
    addAndMakeVisible(dial5);
    addAndMakeVisible(dial6);

    /* Stereo Left */

    /* Left type */
    leftMicType.addItem("Omni-Directional", 1);
    leftMicType.addItem("Cardioid", 2);
    leftMicType.addItem("Figure-8", 3);

    leftMicType.onChange = [this] { leftMicTypeChanged(); };
    leftMicType.setSelectedId(1);
    addAndMakeVisible(leftLabel);
    addAndMakeVisible(leftMicType);


    /* X */
    dial7.setSliderStyle(Slider::SliderStyle::Rotary);
    dial7.setTextBoxStyle(Slider::TextBoxAbove, true, 75, 25);
    dial7.setRange(0, dial1.getValue());
    dial7.setValue(30);
    dial7.setTextValueSuffix("m");
    stereoLeftLocationXYZ.setText("Microphone Stereo Left Location X/Y/Z", dontSendNotification);
    addAndMakeVisible(stereoLeftLocationXYZ);

    /* Y */
    dial8.setSliderStyle(Slider::SliderStyle::Rotary);
    dial8.setTextBoxStyle(Slider::TextBoxAbove, true, 75, 25);
    dial8.setRange(0, dial2.getValue());
    dial8.setValue(30);
    dial8.setTextValueSuffix("m");

    /*  Z */
    dial9.setSliderStyle(Slider::SliderStyle::Rotary);
    dial9.setTextBoxStyle(Slider::TextBoxAbove, true, 75, 25);
    dial9.setRange(0, dial3.getValue());
    dial9.setValue(30);
    dial9.setTextValueSuffix("m");


    /* left mic size */
    micStereoLeftSizeDial.setSliderStyle(Slider::SliderStyle::Rotary);
    micStereoLeftSizeDial.setTextBoxStyle(Slider::TextBoxAbove, true, 75, 25);
    micStereoLeftSizeDial.setRange(0, 50);
    micStereoLeftSizeDial.setValue(10);
    micStereoLeftSizeDial.setTextValueSuffix("m");
    stereoLeftOrientationLabel.setText("Stereo Left Mic Radius/Polar Angle/Azuthimal Angle", dontSendNotification);
    addAndMakeVisible(stereoLeftOrientationLabel);

    /* left polar  */
    micStereoLeftPolarDial.setSliderStyle(Slider::SliderStyle::Rotary);
    micStereoLeftPolarDial.setTextBoxStyle(Slider::TextBoxAbove, true, 75, 25);
    micStereoLeftPolarDial.setRange(0, 360);
    micStereoLeftPolarDial.setValue(0);
    micStereoLeftPolarDial.setTextValueSuffix("degrees");

    /* left azi */
    micStereoLeftAzuthimalDial.setSliderStyle(Slider::SliderStyle::Rotary);
    micStereoLeftAzuthimalDial.setTextBoxStyle(Slider::TextBoxAbove, true, 75, 25);
    micStereoLeftAzuthimalDial.setRange(0, 360);
    micStereoLeftAzuthimalDial.setValue(0);
    micStereoLeftAzuthimalDial.setTextValueSuffix("degrees");


    addAndMakeVisible(micStereoLeftSizeDial);
    addAndMakeVisible(micStereoLeftPolarDial);
    addAndMakeVisible(micStereoLeftAzuthimalDial);

    addAndMakeVisible(dial7);
    addAndMakeVisible(dial8);
    addAndMakeVisible(dial9);

    /* Stereo Right */

    rightMicType.addItem("Omni-Directional", 1);
    rightMicType.addItem("Cardioid", 2);
    rightMicType.addItem("Figure-8", 3);

    rightMicType.onChange = [this] { rightMicTypeChanged(); };
    rightMicType.setSelectedId(1);
    addAndMakeVisible(rightLabel);
    addAndMakeVisible(rightMicType);

    /* X */
    dial10.setSliderStyle(Slider::SliderStyle::Rotary);
    dial10.setTextBoxStyle(Slider::TextBoxAbove, true, 75, 25);
    dial10.setRange(0, dial1.getValue());
    dial10.setValue(50);
    dial10.setTextValueSuffix("m");
    stereoRightLocationXYZ.setText("Microphone Stereo Left Location X/Y/Z", dontSendNotification);
    addAndMakeVisible(stereoRightLocationXYZ);

    /* Y */
    dial11.setSliderStyle(Slider::SliderStyle::Rotary);
    dial11.setTextBoxStyle(Slider::TextBoxAbove, true, 75, 25);
    dial11.setRange(0, dial2.getValue());
    dial11.setValue(50);
    dial11.setTextValueSuffix("m");

    /*  Z */
    dial12.setSliderStyle(Slider::SliderStyle::Rotary);
    dial12.setTextBoxStyle(Slider::TextBoxAbove, true, 75, 25);
    dial12.setRange(0, dial3.getValue());
    dial12.setValue(50);
    dial12.setTextValueSuffix("m");

    /* right mic size */
    micStereoRightSizeDial.setSliderStyle(Slider::SliderStyle::Rotary);
    micStereoRightSizeDial.setTextBoxStyle(Slider::TextBoxAbove, true, 75, 25);
    micStereoRightSizeDial.setRange(0, 50);
    micStereoRightSizeDial.setValue(10);
    micStereoRightSizeDial.setTextValueSuffix("m");
    stereoRightOrientationLabel.setText("Stereo Right Mic Radius/Polar Angle/Azuthimal Angle", dontSendNotification);
    addAndMakeVisible(stereoRightOrientationLabel);

    /* right polar  */
    micStereoRightPolarDial.setSliderStyle(Slider::SliderStyle::Rotary);
    micStereoRightPolarDial.setTextBoxStyle(Slider::TextBoxAbove, true, 75, 25);
    micStereoRightPolarDial.setRange(0, 360);
    micStereoRightPolarDial.setValue(0);
    micStereoRightPolarDial.setTextValueSuffix("degrees");

    /* right azi */
    micStereoRightAzuthimalDial.setSliderStyle(Slider::SliderStyle::Rotary);
    micStereoRightAzuthimalDial.setTextBoxStyle(Slider::TextBoxAbove, true, 75, 25);
    micStereoRightAzuthimalDial.setRange(0, 360);
    micStereoRightAzuthimalDial.setValue(0);
    micStereoRightAzuthimalDial.setTextValueSuffix("degrees");


    addAndMakeVisible(micStereoRightSizeDial);
    addAndMakeVisible(micStereoRightPolarDial);
    addAndMakeVisible(micStereoRightAzuthimalDial);

    addAndMakeVisible(dial10);
    addAndMakeVisible(dial11);
    addAndMakeVisible(dial12);



    /* Additional settings, samples per ray, max ray*/
    /* ThreeDimensionalRoom room, int numRaysPerSample, int maxRaysInExistence, int numDirectionsPerRay */
    /* 3 additonal dials. */

    numSamplesPerRayLabel.setText("Number of samples per set of Rays Fired.", dontSendNotification);
    addAndMakeVisible(numSamplesPerRayLabel);
    int maxSamplesPerRayValue = 10000;
    //numSamplesPerRay.setText("Number of samples per ray (lower more realistic but slower render times).", dontSendNotification);
    numSamplesPerRayDial13.setSliderStyle(Slider::SliderStyle::LinearHorizontal);
    numSamplesPerRayDial13.setTextBoxStyle(Slider::TextBoxLeft, true, 150, 25);
    numSamplesPerRayDial13.setRange(0, maxSamplesPerRayValue);
    numSamplesPerRayDial13.setValue(1000);
    numSamplesPerRayDial13.setSkewFactorFromMidPoint(maxSamplesPerRayValue / 64);
    numSamplesPerRayDial13.addListener(this);
    addAndMakeVisible(numSamplesPerRayDial13);


    maxRaysInExistenceLabel.setText("Max Rays in existence.", dontSendNotification);
    addAndMakeVisible(maxRaysInExistenceLabel);
    int maxRaysInExistence = 100000;
    //numSamplesPerRay.setText("Number of samples per ray (lower more realistic but slower render times).", dontSendNotification);
    maxRaysInExistenceDial14.setSliderStyle(Slider::SliderStyle::LinearHorizontal);
    maxRaysInExistenceDial14.setTextBoxStyle(Slider::TextBoxLeft, true, 150, 25);
    maxRaysInExistenceDial14.setRange(0, maxRaysInExistence);
    maxRaysInExistenceDial14.setValue(7000);
    maxRaysInExistenceDial14.setSkewFactorFromMidPoint(maxRaysInExistence / 128);
    maxRaysInExistenceDial14.addListener(this);
    addAndMakeVisible(maxRaysInExistenceDial14);


    numberOfRaysPerSampleLabel.setText("How many different directions rays are fired in.", dontSendNotification);
    addAndMakeVisible(numberOfRaysPerSampleLabel);
    int numberOfRaysPerSample = 1000;
    //numSamplesPerRay.setText("Number of samples per ray (lower more realistic but slower render times).", dontSendNotification);
    numberOfRaysPerSampleDial15.setSliderStyle(Slider::SliderStyle::LinearHorizontal);
    numberOfRaysPerSampleDial15.setTextBoxStyle(Slider::TextBoxLeft, true, 150, 25);
    numberOfRaysPerSampleDial15.setRange(0, numberOfRaysPerSample);
    numberOfRaysPerSampleDial15.setValue(75);
    numberOfRaysPerSampleDial15.setSkewFactorFromMidPoint(numberOfRaysPerSample / 128);
    numberOfRaysPerSampleDial15.addListener(this);
    addAndMakeVisible(numberOfRaysPerSampleDial15);

    /* Ray Trace Buttons */


    startRayTraceButton.onClick = [this] { startRayTraceButtonClicked(); };
    startRayTraceButton.setColour(juce::TextButton::buttonColourId, juce::Colours::green);
    startRayTraceButton.setEnabled(false);
    addAndMakeVisible(&startRayTraceButton);

    
    

    formatManager.registerBasicFormats();
    transport.addChangeListener(this);

    // Some platforms require permissions to open input channels so request that here
    if (juce::RuntimePermissions::isRequired (juce::RuntimePermissions::recordAudio)
        && ! juce::RuntimePermissions::isGranted (juce::RuntimePermissions::recordAudio))
    {
        juce::RuntimePermissions::request (juce::RuntimePermissions::recordAudio,
                                           [&] (bool granted) { setAudioChannels (granted ? 2 : 0, 2); });
    }
    else
    {
        // Specify the number of input and output channels that we want to open
        setAudioChannels (2, 2);
    }
}



void MainComponent::sliderValueChanged(juce::Slider* slider)
{
    /* Update Sound Source Dials */
    dial4.setRange(0, dial1.getValue());
    dial5.setRange(0, dial2.getValue());
    dial6.setRange(0, dial3.getValue());

    /* Update Stereo Left Dials */
    dial7.setRange(0, dial1.getValue());
    dial8.setRange(0, dial2.getValue());
    dial9.setRange(0, dial3.getValue());
    
    /* Update Stereo Right Dials */
    dial10.setRange(0, dial1.getValue());
    dial11.setRange(0, dial2.getValue());
    dial12.setRange(0, dial3.getValue());

    d1_RoomSizeX = dial1.getValue();
    d2_RoomSizeY = dial2.getValue();
    d3_RoomSizeZ = dial3.getValue();

    d4_SoundSourceLocX = dial4.getValue();
    d5_SoundSourceLocY = dial5.getValue();
    d6_SoundSourceLocZ = dial6.getValue();

    d7_MicStereoLeftLocX = dial7.getValue();
    d8_MicStereoLeftLocY = dial8.getValue();
    d9_MicStereoLeftLocZ = dial9.getValue();

    micStereoLeftSize = micStereoLeftSizeDial.getValue();
    micStereoLeftPolar = micStereoLeftPolarDial.getValue();
    micStereoLeftAzuthimal = micStereoLeftAzuthimalDial.getValue();

    d10_MicStereoRightLocX = dial10.getValue();
    d11_MicStereoRightLocY = dial11.getValue();
    d12_MicStereoRightLocZ = dial12.getValue();

    micStereoRightSize = micStereoRightSizeDial.getValue();
    micStereoRightPolar = micStereoRightPolarDial.getValue();
    micStereoRightAzuthimal = micStereoRightAzuthimalDial.getValue();

    gainDialValue = gainDial.getValue();
    
    numSamplesPerRayDial13Value = numSamplesPerRayDial13.getValue();
    maxRaysInExistenceDial14Value = maxRaysInExistenceDial14.getValue();
    numberOfRaysPerSampleDial15Value = numberOfRaysPerSampleDial15.getValue();
    

    if (slider == &gainDial) {
        /* access all samples in the output file and update the gain values */
        AudioFile<double> reverbedAudioGainControl;

        reverbedAudioGainControl.load(URLofAudioFileForGain);

        int sampleRate = reverbedAudioGainControl.getSampleRate();
        int bitDepth = reverbedAudioGainControl.getBitDepth();
        int numSamples = reverbedAudioGainControl.getNumSamplesPerChannel();
        double lengthInSeconds = reverbedAudioGainControl.getLengthInSeconds();
        int numChannels = reverbedAudioGainControl.getNumChannels();
        bool isMono = reverbedAudioGainControl.isMono();
        bool isStereo = reverbedAudioGainControl.isStereo();

        for (int channel = 0; channel < numChannels; channel++)
        {
            for (int i = 0; i < numSamples; i++)
            {
                //DBG("Before "); DBG(audioFile.samples[channel][i]);
                reverbedAudioGainControl.samples[channel][i] = reverbedAudioGainControl.samples[channel][i] * gainDialValue;
                //DBG("After "); DBG(audioFile.samples[channel][i]);
            }
        }

        DBG("URLofAudioFileAfterGain");
        DBG(URLofAudioFileAfterGain);
        std::string tempOut = URLofAudioFile;
        tempOut.append("   Output.wav");
        reverbedAudioGainControl.save(tempOut);
    }


}

/*void MainComponent::roomSizeMenuChanged()
{

    switch (roomSizeMenu.getSelectedId())
    {
    case 1:    break;
    case 2:      break;

    case 3:    break;
    }
}*/

void MainComponent::openButtonClicked()
{
    //console out clicked
    DBG("clicked");

    juce::FileChooser chooser ("Choose an audio file (wav/aiff)", juce::File::getSpecialLocation(juce::File::userMusicDirectory), "*.wav; *.aiff");
    //if user choses file    
    if (chooser.browseForFileToOpen())
    {
        //juce::File myFile;
        auto myFile = chooser.getResult();
        auto myFileURL = chooser.getURLResult();

        auto url = myFileURL.toString(true);
        URLofAudioFile = url.toStdString();




        DBG("before ", URLofAudioFile);
        DBG(URLofAudioFile);
        char absolutePathURL2[11];
        URLofAudioFile.copy(absolutePathURL2, 12, 0);
        absolutePathURL = absolutePathURL2;
        DBG("absolutePathURL");
        absolutePathURL.pop_back();
        DBG(absolutePathURL);
        URLofAudioFile.erase(0, 12);
        std::regex space("%20");
        URLofAudioFile = std::regex_replace(URLofAudioFile, space, " ");
        DBG(URLofAudioFile);
        
        /* Prepares file to use with AudioFile.h */
        audioFile.load(URLofAudioFile);
        if (audioFile.getNumChannels() != 0) {
            startRayTraceButton.setEnabled(true);
        }
        DBG(audioFile.getNumChannels());
        DBG(audioFile.getNumSamplesPerChannel());

        //Returns pointer
        juce::AudioFormatReader* reader = formatManager.createReaderFor(myFile);
        if (reader != nullptr)
        {
            //the Audio file is a set of data in memory when we bring it into audioformatreadersouce so it's going into memory.
            //tempSource is just a pointer.
            auto duration = (float)reader->lengthInSamples / reader->sampleRate;
            std::unique_ptr<juce::AudioFormatReaderSource> tempSource(new juce::AudioFormatReaderSource(reader, true));

            transport.setSource(tempSource.get());

            playSource.reset(tempSource.release());
            playButton.setEnabled(true);
        }
    }
}

void MainComponent::leftMicTypeChanged() {
    switch (leftMicType.getSelectedId()) {
        case 1: micLeftType = 0;
        case 2: micLeftType = 1;
        case 3: micLeftType = 2;
    }
}
void MainComponent::rightMicTypeChanged() {
    switch (rightMicType.getSelectedId()) {
    case 1: micLeftType = 0;
    case 2: micLeftType = 1;
    case 3: micLeftType = 2;
    }
}

void MainComponent::playButtonClicked()
{
    if ((state == Stopped) || (state == Paused)) {
        transportStateChanged(Starting);
    }
    else if (state == Playing) {
        transportStateChanged(Pausing);
    }
}



void MainComponent::stopButtonClicked()
{
    if (state == Paused) {
        transportStateChanged(Stopped);
    }
    else {
        transportStateChanged(Stopping);
    }
}



void MainComponent::startRayTraceButtonClicked()
{
    //Creates Room Object
    ThreeDimensionalRoom room;

    room.width = 70;
    room.height = 70;
    room.length = 70;
    //room.wallDecayFactor = 0.95;

    //draw coordinates in excel or matlab. Write to csv, then graph.
    soundSource sound;
    sound.x = 10;
    sound.y = 10;
    sound.z = 10;
    sound.reverbAmount = 25;
    

    std::vector<microphone> allMics;
    microphone StereoLeft;
    microphone StereoRight;

    StereoLeft.micLocX = d7_MicStereoLeftLocX; StereoLeft.micLocY = d8_MicStereoLeftLocY; StereoLeft.micLocZ = d9_MicStereoLeftLocZ;
    StereoLeft.size = micStereoLeftSize; StereoLeft.micPolarAngle = micStereoLeftPolar; StereoLeft.micAzimuthalAngle = micStereoLeftAzuthimal;
    StereoRight.micLocX = d10_MicStereoRightLocX; StereoRight.micLocY = d11_MicStereoRightLocY; StereoRight.micLocZ = d12_MicStereoRightLocZ;
    StereoRight.size = micStereoRightSize; StereoRight.micPolarAngle = micStereoRightPolar; StereoRight.micAzimuthalAngle = micStereoRightAzuthimal;

    allMics.push_back(StereoLeft);
    allMics.push_back(StereoRight);
    
    sound.allMic = allMics;

    //Here we read in audio data.
    /* https://github.com/adamstark/AudioFile */
    sound.input = audioFile;
    int sampleRate = audioFile.getSampleRate();
    int bitDepth = audioFile.getBitDepth();
    int numSamples = audioFile.getNumSamplesPerChannel();
    double lengthInSeconds = audioFile.getLengthInSeconds();
    int numChannels = audioFile.getNumChannels();
    bool isMono = audioFile.isMono();
    bool isStereo = audioFile.isStereo();

    // or, just use this quick shortcut to print a summary to the console
    audioFile.printSummary();
    //Copy audio file over to input
    sound.input = audioFile;
    //Access the samples directly and create an output buffer
    for (auto& j : sound.allMic) {
        j.buffer.resize(numChannels);
        j.buffer[0].resize(numSamples);
        j.buffer[1].resize(numSamples);
    }
    /* Does the Ray Trace of sound.input */
    //sound.performRayTrace(room, 1000);


    /* Create the Output Buffer */
    AudioFile<double>::AudioBuffer outputBuffer;
    outputBuffer.resize(numChannels);
    outputBuffer[0].resize(numSamples);
    outputBuffer[1].resize(numSamples);

    AudioFile<double> output;
    output.setBitDepth(bitDepth);
    output.setNumChannels(numChannels);
    output.setNumSamplesPerChannel(numSamples);
    output.setSampleRate(sampleRate);

    float distanceRayTravelsPerIterationToEnsureSpeedOfSound;

    distanceRayTravelsPerIterationToEnsureSpeedOfSound = 343 / (sampleRate / numSamplesPerRayDial13Value);

    DBG(distanceRayTravelsPerIterationToEnsureSpeedOfSound);

    /*    numSamplesPerRayDial13Value = numSamplesPerRayDial13.getValue();
    maxRaysInExistenceDial14Value = maxRaysInExistenceDial14.getValue();
    numberOfRaysPerSampleDial15Value = numberOfRaysPerSampleDial15.getValue();
    
    ThreeDimensionalRoom room, int numRaysPerSample, int maxRaysInExistence, int numDirectionsPerRay, float distanceRayTravelled*/

    /* starts the Ray Trace of sound.input */
    //sound.performRayTrace(room, 1000, 7000, 75, distanceRayTravelsPerIterationToEnsureSpeedOfSound);
    //numberOfRaysPerSampleDial15Value;
    DBG(numberOfRaysPerSampleDial15Value);
    DBG(maxRaysInExistenceDial14Value);
    DBG(numSamplesPerRayDial13Value);
    sound.performRayTrace(room, numSamplesPerRayDial13Value, maxRaysInExistenceDial14Value, numberOfRaysPerSampleDial15Value, distanceRayTravelsPerIterationToEnsureSpeedOfSound);



    /* Reconstruct the two buffers. allMic[0] written to StereoLeft Buffer */
    /* allMic[0] written to StereoLeft Buffer */
    outputBuffer = sound.createStereoAudio(sound.allMic[0].buffer, sound.allMic[1].buffer, outputBuffer);
    

    //ofstream MyFile("left.txt");

    
    

    /* Puts Output into Audio File */
    bool ok = output.setAudioBuffer(outputBuffer);
    


    /* This Is the intermediary wav file, that should be read back into JUCE so the user can listen to the audio in app */
    
    URLofAudioFileJUCE = URLofAudioFile;
    for (int i = 0; i < 4; i++) {
        URLofAudioFileJUCE.pop_back();
    }
    URLofAudioFileJUCE.append(" JUCE.wav");
    output.save(URLofAudioFileJUCE);

    std::regex space(" ");
    URLofAudioFileForGain = URLofAudioFileJUCE;
    URLofAudioFileAfterGain = URLofAudioFileJUCE;
    URLofAudioFileJUCE = std::regex_replace(URLofAudioFileJUCE, space, "%20");
    absolutePathURL.pop_back();
    URLofAudioFileToReadBackIntoJUCE = absolutePathURL.append(URLofAudioFileJUCE);

    DBG(URLofAudioFileToReadBackIntoJUCE);

    //URLofAudioFileToReadBackIntoJUCE.erase(0, 8);

    //URL(URLofAudioFileJUCE);
    DBG(URLofAudioFileToReadBackIntoJUCE);

    //URL url2 = URLofAudioFileToReadBackIntoJUCE;
    URL url = URL(URLofAudioFileToReadBackIntoJUCE);
    File audioFile = url.getLocalFile();

    AudioFormatManager formatManagerOutput;
    formatManagerOutput.registerBasicFormats();


    //Returns pointer
    juce::AudioFormatReader* readerOutput = formatManagerOutput.createReaderFor(audioFile);
    //ScopedPointer<AudioFormatReader> reader = formatManager.createReaderFor(file);
    if (readerOutput != nullptr)
    {
        //the Audio file is a set of data in memory when we bring it into audioformatreadersouce so it's going into memory.
        //tempSource is just a pointer.
        DBG("notNull");
        auto duration = (float)readerOutput->lengthInSamples / readerOutput->sampleRate;
        std::unique_ptr<juce::AudioFormatReaderSource> tempSource(new juce::AudioFormatReaderSource(readerOutput, true));

        transport.setSource(tempSource.get());
        //transportStateChanged(Stopped);

        playSource.reset(tempSource.release());

        //playButtonOutput.setEnabled(true);

    }
    

    //output.save("OutputWav.wav", AudioFileFormat::Wave);

    //Aiff file
    //output.save("OutputAiff.aif", AudioFileFormat::Aiff);

}

void MainComponent::transportStateChanged(TransportState newState)
{
    if (state != newState)
    {
        state = newState;

        switch (state)
        {
        case Stopped:
            playButton.setButtonText("Play");
            stopButton.setButtonText("Stop");
            stopButton.setEnabled(false);
            transport.setPosition(0.0);
            break;

        case Starting:
            transport.start();
            break;

        case Playing:
            playButton.setButtonText("Pause");
            stopButton.setButtonText("Stop");
            stopButton.setEnabled(true);
            break;

        case Pausing:
            transport.stop();
            break;

        case Paused:
            playButton.setButtonText("Resume");
            stopButton.setButtonText("Return to Zero");
            break;

        case Stopping:
            transport.stop();
            break;
        }
    }
}


void MainComponent::changeListenerCallback(ChangeBroadcaster* source)
{
    if (source == &transport)
    {
        if (transport.isPlaying())
            transportStateChanged(Playing);
        else if ((state == Stopping) || (state == Playing))
            transportStateChanged(Stopped);
        else if (Pausing == state)
            transportStateChanged(Paused);
    }
}
/*void MainComponent::startRayTraceButtonClicked()
{
}*/

MainComponent::~MainComponent()
{
    // This shuts down the audio device and clears the audio source.
    shutdownAudio();
}

//==============================================================================
void MainComponent::prepareToPlay (int samplesPerBlockExpected, double sampleRate)
{
    // This function will be called when the audio device is started, or when
    // its settings (i.e. sample rate, block size, etc) are changed.

    // You can use this function to initialise any resources you might need,
    // but be careful - it will be called on the audio thread, not the GUI thread.

    // For more details, see the help for AudioProcessor::prepareToPlay()

    transport.prepareToPlay(samplesPerBlockExpected, sampleRate);
    //transportOutput.prepareToPlay(samplesPerBlockExpected, sampleRate);
}


void MainComponent::getNextAudioBlock (const juce::AudioSourceChannelInfo& bufferToFill)
{
    bufferToFill.clearActiveBufferRegion();
    transport.getNextAudioBlock(bufferToFill);
    bufferToFill.startSample;


    for (auto channel = 0; channel < bufferToFill.buffer->getNumChannels(); ++channel)
    {
        auto* buffer = bufferToFill.buffer->getWritePointer(channel, bufferToFill.startSample);
        for (auto sample = 0; sample < bufferToFill.numSamples; ++sample)
        {
            buffer[sample] = buffer[sample] * gainDialValue;
        }
    }
}

void MainComponent::releaseResources()
{
    // This will be called when the audio device stops, or when it is being
    // restarted due to a setting change.

    // For more details, see the help for AudioProcessor::releaseResources()
}

//==============================================================================
void MainComponent::paint (juce::Graphics& g)
{
    // (Our component is opaque, so we must completely fill the background with a solid colour)
    g.fillAll (getLookAndFeel().findColour (juce::ResizableWindow::backgroundColourId));

    // You can add your drawing code here!
}

void MainComponent::resized()
{
    // This is called when the MainContentComponent is resized.
    // If you add any child components, this is where you should
    // update their positions.

    const int border = 20;

    openButton.setBounds(10, 10, 80, 30);
    playButton.setBounds(110, 10, 80, 30);
    stopButton.setBounds(210, 10, 80, 30);


    const int dialWidth = 85;
    const int dialHeight = 85;

    /* Settings For Ray Tracer */
    playbackSettingsLabel.setBounds(border + 5 * dialWidth, 5, 3 * dialWidth, dialHeight);
    audioSettings->setBounds(210, 50, getWidth() - 210, 30);



    gainLabel.setBounds(border + dialWidth - 20, 25, dialWidth, dialHeight);
    gainDial.setBounds(border, 75, 2*dialWidth, 2*dialHeight);

    roomSizeXYZ.setBounds(border + dialWidth, 240, 3 * dialWidth, dialHeight);
    dial1.setBounds(border+dialWidth, 290, dialWidth, dialHeight);
    dial2.setBounds(border+2*(dialWidth), 290, dialWidth, dialHeight);
    dial3.setBounds(border+3*(dialWidth), 290, dialWidth, dialHeight);

    soundSourceLocationXYZ.setBounds(border + 5*dialWidth, 240, 3 * dialWidth, dialHeight);
    dial4.setBounds(border + 5*dialWidth, 290, dialWidth, dialHeight);
    dial5.setBounds(border + 6*(dialWidth), 290, dialWidth, dialHeight);
    dial6.setBounds(border + 7*(dialWidth), 290, dialWidth, dialHeight);


    leftLabel.setBounds(border + dialWidth, 380, getWidth() - 20, 30);
    leftMicType.setBounds(border + dialWidth, 420, dialWidth * 3, 20);

    //rightLabel.setBounds(border + 5*dialWidth, 410, getWidth() - 20, 20);
    rightMicType.setBounds(border + 5*dialWidth, 420, dialWidth*3, 20);

    stereoLeftLocationXYZ.setBounds(border + dialWidth, 430, 3 * dialWidth, dialHeight);
    dial7.setBounds(border + dialWidth, 480, dialWidth, dialHeight);
    dial8.setBounds(border + 2 * (dialWidth), 480, dialWidth, dialHeight);
    dial9.setBounds(border + 3 * (dialWidth), 480, dialWidth, dialHeight);

    stereoLeftOrientationLabel.setBounds(border + dialWidth, 560, 3 * dialWidth, dialHeight);;
    micStereoLeftSizeDial.setBounds(border + dialWidth, 610, dialWidth, dialHeight);;
    micStereoLeftPolarDial.setBounds(border + 2 * (dialWidth), 610, dialWidth, dialHeight);;
    micStereoLeftAzuthimalDial.setBounds(border + 3 * (dialWidth), 610, dialWidth, dialHeight);;

    stereoRightLocationXYZ.setBounds(border + 5 * dialWidth, 430, 3 * dialWidth, dialHeight);
    dial10.setBounds(border + 5 * dialWidth, 480, dialWidth, dialHeight);
    dial11.setBounds(border + 6 * (dialWidth), 480, dialWidth, dialHeight);
    dial12.setBounds(border + 7 * (dialWidth), 480, dialWidth, dialHeight);

    stereoRightOrientationLabel.setBounds(border + 5 * dialWidth, 560, 3 * dialWidth, dialHeight);;
    micStereoRightSizeDial.setBounds(border + 5 * dialWidth, 610, dialWidth, dialHeight);;
    micStereoRightPolarDial.setBounds(border + 6 * (dialWidth), 610, dialWidth, dialHeight);;
    micStereoRightAzuthimalDial.setBounds(border + 7 * (dialWidth), 610, dialWidth, dialHeight);;

    numSamplesPerRayLabel.setBounds(border + (dialWidth), 685, dialWidth * 7, dialHeight/2);
    numSamplesPerRayDial13.setBounds(border + (dialWidth), 685, dialWidth*7, dialHeight);

    maxRaysInExistenceLabel.setBounds(border + (dialWidth), 730, dialWidth * 7, dialHeight/2);
    maxRaysInExistenceDial14.setBounds(border + (dialWidth), 730, dialWidth*7, dialHeight);

    numberOfRaysPerSampleLabel.setBounds(border + (dialWidth), 775, dialWidth * 7, dialHeight/2);
    numberOfRaysPerSampleDial15.setBounds(border + (dialWidth), 775, dialWidth*7, dialHeight);

    //roomSizeMenu.setBounds(10, 350, getWidth() - 20, 30);


    startRayTraceButton.setBounds(10, 850, getWidth() - 20, 30);
}



