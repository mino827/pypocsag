# pypocsag

POCSAG made easy.

This code was written as a solution to [0xoposec challenge](https://www.meetup.com/0xOPOSEC/events/255388262/).

Raw file can be found [here](https://github.com/zezadas/pypocsag/blob/master/challenge/sdr.raw) and it was captured at a 1Mbps of sample rate and the signal could be found at 153.35 Mhz on the RF spectrum.

My write up can be found [here](https://github.com/zezadas/pypocsag/blob/master/challenge/0x6F_0xOPOSEC_CHALLENGE.pdf).

A python implementation with minimal dependencies.

This pocsag solution was highly based on this specification [http://www.braddye.com/pocsag.html](http://www.braddye.com/pocsag.html).

Thanks to Nuno Humberto for BCH implementation.

## How to use:
Just run it on your terminal using python.
```
python getMessage.py binary_file
```
Input file *binary_file* is a text file containing message bits( 1s and 0s).

## How do I get the message bits from a rawIQ?:
This is not easy to explain by text. Just use inspectrum and a few google skills.

## Why this implementation of pocsag?:
I have written this code has an alternative to decode pocsag messages. This solution has minimal dependencies, oposing to other solutions that demand gqrx,multimon,nc,sox,etc.

## How can I test this decoder?:
Create a message using [gr-mixalot](https://github.com/unsynchronized/gr-mixalot) and gnuradio-companion.
If it fails to compile check my [patch](https://aur.archlinux.org/packages/gr-mixalot-git/).


