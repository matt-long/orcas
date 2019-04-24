import numpy as np
import calendar

month_sname = [calendar.month_name[i+1][0:3] for i in range(12)]
month_letter = [calendar.month_name[i+1][0] for i in range(12)]

dpm = {'noleap': np.array([31.,28.,31.,30.,31.,30.,31.,31.,30.,31.,30.,31.]),
       '365_day': np.array([31.,28.,31.,30.,31.,30.,31.,31.,30.,31.,30.,31.]),
       'standard': np.array([31.,28.,31.,30.,31.,30.,31.,31.,30.,31.,30.,31.]),
       'gregorian': np.array([31.,28.,31.,30.,31.,30.,31.,31.,30.,31.,30.,31.]),
       'proleptic_gregorian': np.array([31.,28.,31.,30.,31.,30.,31.,31.,30.,31.,30.,31]),
       'all_leap': np.array([31.,29.,31.,30.,31.,30.,31.,31.,30.,31.,30.,31.]),
       '366_day': np.array([31.,29.,31.,30.,31.,30.,31.,31.,30.,31.,30.,31.]),
       '360_day': np.array([30.,30.,30.,30.,30.,30.,30.,30.,30.,30.,30.,30.])}

def leap_year(year, calendar='standard'):
    """Determine if year is a leap year"""
    leap = False
    if ((calendar in ['standard', 'gregorian',
        'proleptic_gregorian', 'julian']) and
        (year % 4 == 0)):
        leap = True
        if ((calendar == 'proleptic_gregorian') and
            (year % 100 == 0) and
            (year % 400 != 0)):
            leap = False
        elif ((calendar in ['standard', 'gregorian']) and
                 (year % 100 == 0) and (year % 400 != 0) and
                 (year < 1583)):
            leap = False
    return leap

def get_dpm(time, calendar='standard'):
    """
    return an array of days per month corresponding to the months provided in `months`
    """
    month_length = np.zeros(len(time), dtype=np.float)

    cal_days = dpm[calendar]

    for i, (month, year) in enumerate(zip(time.month, time.year)):
        month_length[i] = cal_days[month-1]
        if leap_year(year, calendar=calendar) and month==2:
            month_length[i] += 1
    return month_length

