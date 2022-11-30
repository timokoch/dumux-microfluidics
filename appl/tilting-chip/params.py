"""
Base parameter sets for different scenarios
These are parameters to be passed to the C++ program as command line arguments
"""

BIG_RESERVOIR = {
    "Grid.File": "big-extruded.msh",
    "Problem.InitialVolumeInMicroLiter": "300",
    "Problem.MeasurementPoint1": "-8.2 -10.5 -0.49",
    "Problem.MeasurementPoint2": "-8.2 10.5 -0.49",
}

SMALL_RESERVOIR = {
    "Grid.File": "small-extruded.msh",
    "Problem.InitialVolumeInMicroLiter": "200",
    "Problem.MeasurementPoint1": "-4.9 -6.5 -0.49",
    "Problem.MeasurementPoint2": "-4.9 6.5 -0.49",
}

HIGH_ASPECT_RATIO_CHANNEL_LONG = {
    "Problem.DryingThresholdInMilliMeter": "0.7",
    "Problem.SingleChannelVolumeInMicroLiter": "26.24",
    "Problem.ChannelLengthInMeter": "16.4e-3",
    "Problem.ChannelAreaInSquareMeter": str(3.2*0.5e-6),
}

HIGH_ASPECT_RATIO_CHANNEL_SHORT = {
    "Problem.DryingThresholdInMilliMeter": "0.7",
    "Problem.SingleChannelVolumeInMicroLiter": "15.68",
    "Problem.ChannelLengthInMeter": "9.8e-3",
    "Problem.ChannelAreaInSquareMeter": str(3.2*0.5e-6),
}

LOW_ASPECT_RATIO_CHANNEL_LONG = {
    "Problem.DryingThresholdInMilliMeter": "0.5",
    "Problem.SingleChannelVolumeInMicroLiter": "15.74",
    "Problem.ChannelLengthInMeter": "16.4e-3",
    "Problem.ChannelAreaInSquareMeter": str(1.2*0.8e-6),
}

LOW_ASPECT_RATIO_CHANNEL_SHORT = {
    "Problem.DryingThresholdInMilliMeter": "0.5",
    "Problem.SingleChannelVolumeInMicroLiter": "9.41",
    "Problem.ChannelLengthInMeter": "9.8e-3",
    "Problem.ChannelAreaInSquareMeter": str(1.2*0.8e-6),
}

BIG_RESERVOIR_HIGH_ASPECT_CHANNEL = {
    **BIG_RESERVOIR,
    **HIGH_ASPECT_RATIO_CHANNEL_LONG,
}

BIG_RESERVOIR_LOW_ASPECT_CHANNEL = {
    **BIG_RESERVOIR,
    **LOW_ASPECT_RATIO_CHANNEL_LONG,
}

SMALL_RESERVOIR_HIGH_ASPECT_CHANNEL = {
    **SMALL_RESERVOIR,
    **HIGH_ASPECT_RATIO_CHANNEL_SHORT,
}

SMALL_RESERVOIR_LOW_ASPECT_CHANNEL = {
    **SMALL_RESERVOIR,
    **LOW_ASPECT_RATIO_CHANNEL_SHORT,
}

TRANSMISSIBILITY_20_DEGREE = {
    "HIGH_ASPECT_RATIO_CHANNEL_SHORT": {
        "Problem.ChannelTransmissibility": "3.0664e-9",
        "Problem.ChannelWSSFactor": "6.17e+06", # WSS/Q
    },
    "HIGH_ASPECT_RATIO_CHANNEL_LONG": {
        "Problem.ChannelTransmissibility": "1.8323e-9",
        "Problem.ChannelWSSFactor": "6.17e+06", # WSS/Q
    },
    "LOW_ASPECT_RATIO_CHANNEL_SHORT": {
        "Problem.ChannelTransmissibility": "3.0683e-9",
        "Problem.ChannelWSSFactor": "9.67e+06", # WSS/Q
    },
    "LOW_ASPECT_RATIO_CHANNEL_LONG": {
        "Problem.ChannelTransmissibility": "1.8334e-9",
        "Problem.ChannelWSSFactor": "9.67e+06", # WSS/Q
    },
}

TRANSMISSIBILITY_37_DEGREE = {
    "HIGH_ASPECT_RATIO_CHANNEL_SHORT": {
        "Problem.ChannelTransmissibility": "3.8330e-9",
        "Problem.ChannelWSSFactor": "4.93e+6", # WSS/Q
    },
    "HIGH_ASPECT_RATIO_CHANNEL_LONG": {
        "Problem.ChannelTransmissibility": "2.2905e-9",
        "Problem.ChannelWSSFactor": "4.93e+6", # WSS/Q
    },
    "LOW_ASPECT_RATIO_CHANNEL_SHORT": {
        "Problem.ChannelTransmissibility": "3.8353e-9",
        "Problem.ChannelWSSFactor": "7.74e+6", # WSS/Q
    },
    "LOW_ASPECT_RATIO_CHANNEL_LONG": {
        "Problem.ChannelTransmissibility": "2.2918e-9",
        "Problem.ChannelWSSFactor": "7.74e+6", # WSS/Q
    },
}

DEFAULT_KINETIC = {
    "Problem.RotationsPerMinute": "3",
    "Problem.Angle": "19.0",
    "TimeLoop.Cycles": "1",
    "TimeLoop.Dt": "0.1",
}

DEFAULT_FLUX_LIMITER = {
    "Problem.CapillaryStopPressure": "20", # Pa
    "Problem.MaxFluxChangePerSecond": "500", # µl/s/s
}

# Validation with experiments:
# For that we use 20 degrees and the big channel with high aspect ratio
DEFAULT_CASE_VALIDATION = {
    **BIG_RESERVOIR_HIGH_ASPECT_CHANNEL,
    **(TRANSMISSIBILITY_20_DEGREE["HIGH_ASPECT_RATIO_CHANNEL_LONG"]),
    **DEFAULT_KINETIC,
    **DEFAULT_FLUX_LIMITER,
}

DEFAULT_CASE_VISUALIZATION_BIG = DEFAULT_CASE_VALIDATION
DEFAULT_CASE_VISUALIZATION_SMALL = {
    **SMALL_RESERVOIR_HIGH_ASPECT_CHANNEL,
    **(TRANSMISSIBILITY_20_DEGREE["HIGH_ASPECT_RATIO_CHANNEL_SHORT"]),
    **DEFAULT_KINETIC,
    **DEFAULT_FLUX_LIMITER,
}

# WSS and Q calculations:
# For that we use 37 degrees and both channels with low aspect ratio
DEFAULT_CASE_WSS_BIG = {
    **BIG_RESERVOIR_LOW_ASPECT_CHANNEL,
    **(TRANSMISSIBILITY_37_DEGREE["LOW_ASPECT_RATIO_CHANNEL_LONG"]),
    **DEFAULT_KINETIC,
    **DEFAULT_FLUX_LIMITER,
}

DEFAULT_CASE_WSS_SMALL = {
    **SMALL_RESERVOIR_LOW_ASPECT_CHANNEL,
    **(TRANSMISSIBILITY_37_DEGREE["LOW_ASPECT_RATIO_CHANNEL_SHORT"]),
    **DEFAULT_KINETIC,
    **DEFAULT_FLUX_LIMITER,
}


def make_cmd(args, output_filename, prefix_path=None):
    """Create command with Dumux command line arguments"""
    cmd = ["./flow_model",] if prefix_path is None else [f"{prefix_path}/flow_model",]
    for k, v in args.items():
        cmd += [f"-{k}", f"{v}",]
    cmd += ["-Problem.OutputFileName", output_filename,]
    return cmd
