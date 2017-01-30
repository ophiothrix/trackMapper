#!/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$DIR"

source MapNormTracks.param.txt

R --vanilla --args $genes $window $upperlim $upstream $downstream $pathToBam $useNormalisedLibrarySize < MapNormTracksNoGrouping.R 
