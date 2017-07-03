#!/bin/bash
while [ ! -f .executed ]; do
  run focal_point
  sleep 10
done
