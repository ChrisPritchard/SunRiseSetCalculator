# SunRiseSetCalculator
A small set of classes which take a latitude/longitude, date and a timezone, and calculate the approximate times for sun rise and set.

The timezone requirement might be something you think can be auto-calculated based on lat long, but you'd be wrong. Or at least its much more difficult than the code for calculating the sun times, largely because timezones are to a degree arbitary (they are not in fixed-size bands across the world, but wriggle a bit based on political boundaries).

Accuracy has been tested in the NZ timezone (roughly GMT+12), and despite being a calculation based largely on the worlds geodesic projection relative to the sun, tends to be withing 5 or 10 minutes of other estimates I have seen.

While the code is packaged as NetStandard 1.4, the classes should be fully compatible with .NET Framework 4.6 (or earlier, if you convert some of the C# 6 statements to their older equivalents).