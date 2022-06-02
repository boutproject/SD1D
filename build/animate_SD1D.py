from xbout import open_boutdataset
import matplotlib.pyplot as plt
import click


@click.command()
@click.option("-d", default="", help="Directory of BOUT++ output files.")
def animate_SD1D_output(d):
    bd = open_boutdataset(
        f"./{d}/BOUT.dmp.*.nc", inputfilepath=f"./{d}/BOUT.inp"
    ).squeeze(drop=True)
    try:
        bd.bout.animate_list(["Ne", "NVi", "P", "Nn", "NVn", "Pn"], nrows=2, ncols=3)
    except Exception:
        plt.close()
        bd.bout.animate_list(["Ne", "NVi", "P"], nrows=3, ncols=1)
    plt.show()


if __name__ == "__main__":
    animate_SD1D_output()
