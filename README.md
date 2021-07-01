This is a set of tools that are being used in the TDCOSMO collaboration to
study the environments around strong gravitational lenses.


## Installation

lenskappa has been tested and is working with python >= 3.7. It may work on older versions of python, but I recommend just creating a new conda or virtualenv environment.

Once this environment has been created, you can install lenskappa by doing the following.

```
git clone https://github.com/PatrickRWells/lenskappa.git  
cd lenskappa  
pip install -r requirements.txt  
pip install -e .  
```
Since lenskappa is in early development, it is recommended you install it in development mode (denoted by the `-e` flag in the install statement). You can then update the package simply by pulling the latest changes from github.

## Usage

### Adding control field data

Lenskappa can keep track of control field data for you, but you first have to add it using the lenskappa_add_surveydata script, which is installed automatically along with the package.

As an example, suppose you were interested in computing weighted number counts for the W02 field of the HSC survey. To add the W02 catalog to lenskappa, you would do the following.

```
lenskappa_add_surveydata hsc W02 catalog /path/to/your/data/W02.csv
```
### Computing weighted number counts

Once you have added your control field data, you can compute weight ratios for any system of interest. See example.py for an example of this.
