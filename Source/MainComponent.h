#pragma once
#include <JuceHeader.h>
#include "AudioFile.h"
using namespace juce;



//==============================================================================
/*
    This component lives inside our window, and this is where you should put all
    your controls and content.
*/
class MainComponent  : public AudioAppComponent,
                       public ChangeListener,
                       public Slider::Listener
{
public:
    //==============================================================================
    MainComponent();
    ~MainComponent() override;

    //==============================================================================
    void prepareToPlay (int samplesPerBlockExpected, double sampleRate) override;
    void getNextAudioBlock (const AudioSourceChannelInfo& bufferToFill) override;
    void releaseResources() override;

    //==============================================================================
    void paint (Graphics& g) override;
    void resized() override;

private:
    AudioDeviceManager otherDeviceManager;
    std::unique_ptr<AudioDeviceSelectorComponent> audioSettings;
    

    enum TransportState
    {
        /* Open File States*/
        Stopped,
        Starting,
        Stopping,
        Playing,
        Pausing,
        Paused,

    };

    TransportState state;

    void openButtonClicked();

    void playButtonClicked();
    void stopButtonClicked();

    void leftMicTypeChanged();
    void rightMicTypeChanged();

    void startRayTraceButtonClicked();

    void transportStateChanged(TransportState newState);

    //virtual function
    void changeListenerCallback(ChangeBroadcaster* source) override;

    AudioFormatManager formatManager;
    AudioFormatManager formatManagerOutput;

    std::unique_ptr<AudioFormatReaderSource> playSource;
    AudioTransportSource transport;


    AudioFile<double> audioFile;

    std::string URLofAudioFile;
    std::string absolutePathURL;

    std::string URLofAudioFileForGain;
    std::string URLofAudioFileAfterGain;

    std::string URLofAudioFileJUCE;
    std::string URLofAudioFileToReadBackIntoJUCE;


    int d1_RoomSizeX;
    int d2_RoomSizeY;
    int d3_RoomSizeZ;

    int d4_SoundSourceLocX;
    int d5_SoundSourceLocY;
    int d6_SoundSourceLocZ;
    
    int micLeftType;
    int d7_MicStereoLeftLocX;
    int d8_MicStereoLeftLocY;
    int d9_MicStereoLeftLocZ;

    int micStereoLeftSize;
    int micStereoLeftPolar;
    int micStereoLeftAzuthimal;

    int micRightType;
    int d10_MicStereoRightLocX;
    int d11_MicStereoRightLocY;
    int d12_MicStereoRightLocZ;

    int micStereoRightSize;
    int micStereoRightPolar;
    int micStereoRightAzuthimal;

    int gainSlider;
    float gainDialValue;

    float numSamplesPerRayDial13Value;
    float maxRaysInExistenceDial14Value;
    float numberOfRaysPerSampleDial15Value;


    int position;


    TextButton openButton;
    TextButton playButton;
    TextButton stopButton;


    /* Set up the ray tracer settings */
    juce::Font textFont{ 12.0f };

    //juce::ComboBox roomSizeMenu;
    juce::ComboBox volumeOfRays;

    TextButton startRayTraceButton;
    TextButton playButtonOutput;
    TextButton stopButtonOutput;

    Label playbackSettingsLabel;

    Label gainLabel;
    Slider gainDial;
    int maxGainDialSize;


    int maxRoomSize;
    Label roomSizeXYZ;
    Slider dial1;
    Slider dial2;
    Slider dial3;

    Label soundSourceLocationXYZ;
    Slider dial4;
    Slider dial5;
    Slider dial6;


    Label leftLabel{ {}, "Stereo Left and Right(respectively) Pickup Patterns. (Omni-directional is recommended) \n Polar & Azuthimal are Spherical Coordinates." };
    ComboBox leftMicType;
    
    Label rightLabel{ {}, "Stereo Right Pickup Pattern. (Omni-directional is recommended)" };
    ComboBox rightMicType;

    Label stereoLeftLocationXYZ;
    Slider dial7;
    Slider dial8;
    Slider dial9;

    Label stereoLeftOrientationLabel;
    Slider micStereoLeftSizeDial;
    Slider micStereoLeftPolarDial;
    Slider micStereoLeftAzuthimalDial;

    Label stereoRightLocationXYZ;
    Slider dial10;
    Slider dial11;
    Slider dial12;

    Label stereoRightOrientationLabel;
    Slider micStereoRightSizeDial;
    Slider micStereoRightPolarDial;
    Slider micStereoRightAzuthimalDial;

    Label numSamplesPerRayLabel;
    Slider numSamplesPerRayDial13;

    Label maxRaysInExistenceLabel;
    Slider maxRaysInExistenceDial14;

    Label numberOfRaysPerSampleLabel;
    Slider numberOfRaysPerSampleDial15;

    void sliderValueChanged(juce::Slider* slider);

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (MainComponent)
};
