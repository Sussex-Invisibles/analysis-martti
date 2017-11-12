#!/bin/bash
while [ ! -f .executed ]; do
  run fibre_validation
  sleep 10
done
