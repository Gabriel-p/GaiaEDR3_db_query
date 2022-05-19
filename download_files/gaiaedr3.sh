#!/bin/bash


# Session Name
session="gaiaEDR3"

# Start New Session with our name
tmux new-session -d -s $session

tmux send-keys "exec ./gaia_base.sh" 'C-m'
