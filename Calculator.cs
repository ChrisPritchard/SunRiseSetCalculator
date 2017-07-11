using System;

namespace SunRiseSetCalculator
{
    public class Calculator
    {
        public static TimeSet FindSunRiseAndSet(DateTime dateTime, double latitude, double longitude, TimeZoneInfo timeZoneInfo, bool forceRecalculate = false)
        {
            var utcOffset = timeZoneInfo.GetUtcOffset(dateTime).TotalHours;
            return CalculateSunRiseAndSet(dateTime, utcOffset, longitude, latitude);
        }

        private static TimeSet CalculateSunRiseAndSet(DateTime dt, double tz, double glong, double glat)
        {
            double mjd = GetJulianDate(dt, 0);
            double sglat, cglat, date, ym, yz, utrise = 0, utset = 0;
            double yp, nz, hour, xe, ye, z1, z2, rads = 0.0174532925;

            bool rise, sett, above;
            double[] quadout = new double[4];

            var sinho = Math.Sin(rads * -0.833);
            sglat = Math.Sin(rads * glat);
            cglat = Math.Cos(rads * glat);
            date = mjd - tz / 24;

            rise = false;
            sett = false;
            above = false;
            hour = 1.0;

            ym = SineOfAltitude(2, date, hour - 1.0, glong, cglat, sglat) - sinho;
            if (ym > 0.0) above = true;
            //
            // the while loop finds the sin(alt) for sets of three consecutive
            // hours, and then tests for a single zero crossing in the interval
            // or for two zero crossings in an interval or for a grazing event
            // The flags rise and sett are set accordingly
            //
            nz = 0; z1 = 0; z2 = 0; xe = 0; ye = 0;
            while (hour < 25 && (sett == false || rise == false))
            {
                yz = SineOfAltitude(2, date, hour, glong, cglat, sglat) - sinho;
                yp = SineOfAltitude(2, date, hour + 1.0, glong, cglat, sglat) - sinho;
                quadout = Quad(ym, yz, yp);
                nz = quadout[0];
                z1 = quadout[1];
                z2 = quadout[2];
                xe = quadout[3];
                ye = quadout[4];

                // case when one event is found in the interval
                if (nz == 1)
                {
                    if (ym < 0.0)
                    {
                        utrise = hour + z1;
                        rise = true;
                    }
                    else
                    {
                        utset = hour + z1;
                        sett = true;
                    }
                } // end of nz = 1 case

                // case where two events are found in this interval
                // (rare but whole reason we are not using simple iteration)
                if (nz == 2)
                {
                    if (ye < 0.0)
                    {
                        utrise = hour + z2;
                        utset = hour + z1;
                    }
                    else
                    {
                        utrise = hour + z1;
                        utset = hour + z2;
                    }
                } // end of nz = 2 case

                // set up the next search interval
                ym = yp;
                hour += 2.0;

            } // end of while loop
            //
            // now search has completed, we compile the string to pass back
            // to the main loop. The string depends on several combinations
            // of the above flag (always above or always below) and the rise
            // and sett flags
            //

            var sunRiseDate = ((rise == true || sett == true)) ? (rise == true) ? ConvertToDateTime(dt, utrise) : DateTime.MinValue : (above == true) ? DateTime.MaxValue : DateTime.MaxValue.AddDays(-1);
            var sunSetDate = ((rise == true || sett == true)) ? (sett == true) ? ConvertToDateTime(dt, utset) : DateTime.MinValue.AddDays(1) : (above == true) ? DateTime.MaxValue : DateTime.MaxValue.AddDays(-1);
            return new TimeSet
            {
                SunRise = sunRiseDate.TimeOfDay,
                SunSet = sunSetDate.TimeOfDay
            };
        }

        /// <summary>
        /// returns the integer part - like int() in basic
        /// </summary>
        private static int IntegerPart(double x)
        {
            double a;
            if (x > 0)
            {
                a = Math.Floor(x);
            }
            else
            {
                a = Math.Ceiling(x);
            }
            return Convert.ToInt32(a);
        }

        /// <summary>
        /// returns the fractional part of x as used in minimoon and minisun
        /// </summary>
        private static double Frac(double x)
        {
            double a;
            a = x - Math.Floor(x);
            if (a < 0) a += 1;
            return a;
        }

        /// <summary>
        /// returns an angle in degrees in the range 0 to 360
        /// </summary>
        private static double Range(double x)
        {
            double a;
            double b;
            b = x / 360;
            a = 360 * (b - IntegerPart(b));
            if (a < 0)
            {
                a = a + 360;
            }
            return a;
        }

        /// <summary>
        /// Takes the day, month, year and hours in the day and returns the
        /// modified julian day number defined as mjd = jd - 2400000.5
        /// checked OK for Greg era dates - 26th Dec 02
        /// </summary>
        private static double GetJulianDate(DateTime dt, int hour)
        {
            int Day, Month, Year;

            Day = dt.Day;
            Month = dt.Month;
            Year = dt.Year;

            double a, b;
            if (Month <= 2)
            {
                Month = Month + 12;
                Year = Year - 1;
            }
            a = 10000.0 * Year + 100.0 * Month + Day;
            if (a <= 15821004.1)
            {
                b = -2 * Math.Floor((double)(Year + 4716) / 4) - 1179;
            }
            else
            {
                b = Math.Floor((double)Year / 400) - Math.Floor((double)Year / 100) + Math.Floor((double)Year / 4);
            }
            a = 365.0 * Year - 679004.0;
            return (a + b + Math.Floor(30.6001 * (Month + 1)) + Day + hour / 24.0);
        }

        /// <summary>
        /// finds the parabola throuh the three points (-1,ym), (0,yz), (1, yp)
        /// and returns the coordinates of the max/min (if any) xe, ye
        /// the values of x where the parabola crosses zero (roots of the quadratic)
        /// and the number of roots (0, 1 or 2) within the interval [-1, 1]
        ///
        /// well, this routine is producing sensible answers
        ///
        /// results passed as array [nz, z1, z2, xe, ye]
        /// </summary>
        private static double[] Quad(double ym, double yz, double yp)
        {
            double nz, a, b, c, dis, dx, xe, ye, z1 = 0, z2 = 0;
            double[] quadout = new double[5];

            nz = 0;
            a = 0.5 * (ym + yp) - yz;
            b = 0.5 * (yp - ym);
            c = yz;
            xe = -b / (2 * a);
            ye = (a * xe + b) * xe + c;
            dis = b * b - 4.0 * a * c;
            if (dis > 0)
            {
                dx = 0.5 * Math.Sqrt(dis) / Math.Abs(a);
                z1 = xe - dx;
                z2 = xe + dx;
                if (Math.Abs(z1) <= 1.0) nz += 1;
                if (Math.Abs(z2) <= 1.0) nz += 1;
                if (z1 < -1.0) z1 = z2;
            }
            quadout[0] = nz;
            quadout[1] = z1;
            quadout[2] = z2;
            quadout[3] = xe;
            quadout[4] = ye;
            return quadout;
        }

        /// <summary>
        /// Takes the mjd and the longitude (west negative) and then returns
        /// the local sidereal time in hours. Im using Meeus formula 11.4
        /// instead of messing about with UTo and so on
        /// </summary>
        private static double Sidereal(double mjd, double glong)
        {
            double lst, t, d;
            d = mjd - 51544.5;
            t = d / 36525.0;
            lst = Range(280.46061837 + 360.98564736629 * d + 0.000387933 * t * t - t * t * t / 38710000);
            return (lst / 15.0 + glong / 15);
        }

        /// <summary>
        /// returns the ra and dec of the Sun in an array called suneq[]
        /// in decimal hours, degs referred to the equinox of date and using
        /// obliquity of the ecliptic at J2000.0 (small error for +- 100 yrs)
        /// takes t centuries since J2000.0. Claimed good to 1 arcmin
        /// </summary>
        private static double[] MiniSun(double t)
        {
            double p2 = 6.283185307, coseps = 0.91748, sineps = 0.39778;
            double L, M, DL, SL, X, Y, Z, RHO, ra, dec;
            double[] suneq = new double[2]; ;

            M = p2 * Frac(0.993133 + 99.997361 * t);
            DL = 6893.0 * Math.Sin(M) + 72.0 * Math.Sin(2 * M);
            L = p2 * Frac(0.7859453 + M / p2 + (6191.2 * t + DL) / 1296000);
            SL = Math.Sin(L);
            X = Math.Cos(L);
            Y = coseps * SL;
            Z = sineps * SL;
            RHO = Math.Sqrt(1 - Z * Z);
            dec = (360.0 / p2) * Math.Atan(Z / RHO);
            ra = (48.0 / p2) * Math.Atan(Y / (X + RHO));
            if (ra < 0) ra += 24;
            suneq[0] = dec;
            suneq[1] = ra;
            return suneq;
        }

        /// <summary>
        /// takes t and returns the geocentric ra and dec in an array mooneq
        /// claimed good to 5' (angle) in ra and 1' in dec
        /// tallies with another approximate method and with ICE for a couple of dates
        /// </summary>
        private static double[] MiniMoon(double t)
        {
            double p2 = 6.283185307, arc = 206264.8062, coseps = 0.91748, sineps = 0.39778;
            double L0, L, LS, F, D, H, S, N, DL, CB, L_moon, B_moon, V, W, X, Y, Z, RHO, dec, ra;
            double[] mooneq = new double[2];

            L0 = Frac(0.606433 + 1336.855225 * t);  // mean longitude of moon
            L = p2 * Frac(0.374897 + 1325.552410 * t); //mean anomaly of Moon
            LS = p2 * Frac(0.993133 + 99.997361 * t); //mean anomaly of Sun
            D = p2 * Frac(0.827361 + 1236.853086 * t); //difference in longitude of moon and sun
            F = p2 * Frac(0.259086 + 1342.227825 * t); //mean argument of latitude

            // corrections to mean longitude in arcsec
            DL = 22640 * Math.Sin(L);
            DL += -4586 * Math.Sin(L - 2 * D);
            DL += +2370 * Math.Sin(2 * D);
            DL += +769 * Math.Sin(2 * L);
            DL += -668 * Math.Sin(LS);
            DL += -412 * Math.Sin(2 * F);
            DL += -212 * Math.Sin(2 * L - 2 * D);
            DL += -206 * Math.Sin(L + LS - 2 * D);
            DL += +192 * Math.Sin(L + 2 * D);
            DL += -165 * Math.Sin(LS - 2 * D);
            DL += -125 * Math.Sin(D);
            DL += -110 * Math.Sin(L + LS);
            DL += +148 * Math.Sin(L - LS);
            DL += -55 * Math.Sin(2 * F - 2 * D);

            // simplified form of the latitude terms
            S = F + (DL + 412 * Math.Sin(2 * F) + 541 * Math.Sin(LS)) / arc;
            H = F - 2 * D;
            N = -526 * Math.Sin(H);
            N += +44 * Math.Sin(L + H);
            N += -31 * Math.Sin(-L + H);
            N += -23 * Math.Sin(LS + H);
            N += +11 * Math.Sin(-LS + H);
            N += -25 * Math.Sin(-2 * L + F);
            N += +21 * Math.Sin(-L + F);

            // ecliptic long and lat of Moon in rads
            L_moon = p2 * Frac(L0 + DL / 1296000);
            B_moon = (18520.0 * Math.Sin(S) + N) / arc;

            // equatorial coord conversion - note fixed obliquity
            CB = Math.Cos(B_moon);
            X = CB * Math.Cos(L_moon);
            V = CB * Math.Sin(L_moon);
            W = Math.Sin(B_moon);
            Y = coseps * V - sineps * W;
            Z = sineps * V + coseps * W;
            RHO = Math.Sqrt(1.0 - Z * Z);
            dec = (360.0 / p2) * Math.Atan(Z / RHO);
            ra = (48.0 / p2) * Math.Atan(Y / (X + RHO));
            if (ra < 0) ra += 24;
            mooneq[0] = dec;
            mooneq[1] = ra;
            return mooneq;
        }

        /// <summary>
        /// this rather mickey mouse function takes a lot of
        /// arguments and then returns the sine of the altitude of
        /// the object labelled by iobj. iobj = 1 is moon, iobj = 2 is sun
        /// </summary>
        private static double SineOfAltitude(double iobj, double mjd0, double hour, double glong, double cglat, double sglat)
        {
            double mjd, t, ra, dec, tau, salt, rads = 0.0174532925;
            double[] objpos = new double[2];
            mjd = mjd0 + hour / 24.0;
            t = (mjd - 51544.5) / 36525.0;
            if (iobj == 1)
            {
                objpos = MiniMoon(t);
            }
            else
            {
                objpos = MiniSun(t);
            }
            ra = objpos[1];
            dec = objpos[0];
            // hour angle of object
            tau = 15.0 * (Sidereal(mjd, glong) - ra);
            // sin(alt) of object using the conversion formulas
            salt = sglat * Math.Sin(rads * dec) + cglat * Math.Cos(rads * dec) * Math.Cos(rads * tau);
            return salt;
        }

        private static DateTime ConvertToDateTime(DateTime yearMonthDay, double totalHours)
        {
            var hours = (int)Math.Floor(totalHours);
            var remaining = (totalHours - hours) * 60;
            var minutes = (int)Math.Floor(remaining);
            remaining = (remaining - minutes) * 60;
            var seconds = (int)Math.Floor(remaining);

            return new DateTime(yearMonthDay.Year, yearMonthDay.Month, yearMonthDay.Day, hours, minutes, seconds);
        }
    }
}
