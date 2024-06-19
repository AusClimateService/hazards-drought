import sys

submodule_path = '/g/data/mn51/users/jb6465/drought-github/submodules/gwls'
sys.path.append(submodule_path)

from gwl import get_GWL_syear_eyear


def main():
    result = get_GWL_syear_eyear('CMIP6','ACCESS-ESM1-5','r6i1p1f1','ssp370','1.2')
    print(result)

if __name__ == "__main__":
    main()