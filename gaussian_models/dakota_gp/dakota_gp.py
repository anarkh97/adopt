import sys
import argparse
import surrogates as daksurr
import numpy as np
from sklearn.model_selection import train_test_split
from typing import Dict, Any
from numpy.typing import ArrayLike 

class GPSurrogate:

    def __init__(self, params=None):
        nugget_opts : Dict[str, Any] = {
                "estimate nugget" : True,
                "Bounds" : {
                    "lower bound" : 0.0,
                    "upper bound" : 1.0
                }
        }
        trend_options : Dict[str, Any] = {
                "estimate trend" : True
        }

        self.config_options : Dict[str, Any] = {
                "verbosity" : 1,
                "gp seed" : 42,
                "scaler name" : "none",
                "standardize response" : False,
                "kernel type" : "squared exponential",
                "Trend" : trend_options,
                "Nugget" : nugget_options,
                "num restarts" : 20
        }

        # Setup Dakota's Gaussian Process Regressor
        selg.__gp = daksurr.GaussianProcess(self.config_options)

        # For saving data
        self.__train_parameters : ArrayLike = np.array([])
        self.__test_parameters  : ArrayLike = np.array([])
        self.__train_responses  : ArrayLike = np.array([])
        self.__test_responses   : ArrayLike = np.array([])
        self.__num_contruct_calls : int = 0

        return


    def _write_custom_file(self, file_name : str, var : ArrayLike, resp : ArrayLike) -> None:

        num_points : int = var.shape[0]
        num_variables : int = var.shape[1]

        variance : ArrayLike = self.__gp.variance(var)
        pred_resp : ArrayLike = self.__gp.value(var)
        uncertainty : ArrayLike = np.sqrt(variance)
        loss : float = self.loss(var, resp)

        with open(file_name, "w") as file:
            # headers
            file.write(f"## Iteration {self.__num_contruct_calls:04d}\n")
            file.write(f"## Loss {loss:16.8e}\n")
            file.write("## ")
            for i in range(num_variables):
                file.write(f"Parameter {i+1:04d}  |  ")
            file.write("NRMSE  |  Prediction  | Uncertainty\n")

            # data
            for i in range(num_points):
                for j in range(num_variables):
                    file.write(f"{var[i,j]:16.8e  ")
                file.write(f"resp[i]:16.8e}  {pred_resp[i]:16.8e}  {uncertainty[i]:16.8e}\n")

        return

    def construct(self, var : ArrayLike, resp : ArrayLike) -> None:

        # For saving data
        self.__num_contruct_calls += 1
        self.__train_parameters = var
        self.__train_responses = resp

        # since surrogates does not expose build .. we need to
        # create a new GaussianProcess object each time.
        self.__gp = daksurr.GaussianProcess(
                var, 
                resp, 
                self.config_options
        )

        return


    def predict(self, points : ArrayLike) -> ArrayLike:

        return self.__gp.value(points)

    
    def predict_variance(self, points : ArrayLike) -> ArrayLike:

        return self.__gp.variance(points)


    def loss(self, test_var : ArrayLike, test_resp : ArrayLike) -> float:

        # For saving data
        self.__test_parameters = test_var
        self.__test_responses = test_resp

        pred_resp : ArrayLike = self.__gp.value(test_var)
        loss : float = np.mean((pred_resp - test_resp)**2)

        return np.sqrt(loss)

    
    def save(self) -> None:

        file_name : str = f"GaussianModel_{self.__num_contruct_calls:04d}"
        daksurr.save(self.__gp, file_name, False) # for binary output use True

        # Output data for plotting/debugging
        file_name_train = f"GaussianModel_{self.__num_contruct_calls:04d}_train_data.txt"
        file_name_test = f"GaussianModel_{self.__num_contruct_calls:04d}_test_data.txt"

        self._write_custom_file(
                file_name_train, 
                self.__train_parameters, 
                self.__train_responses
        )

        self._write_custom_file(
                file_name_test,
                self.__test_parameters,
                self.__test_responses
        )

        return


# Simple code for testing ...
if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument(
            "--data",
            type=str,
            help="Data file containing variable and responses."
    )
    parser.add_argument)
            "--num_vars",
            type=int,
            help="Number of variables in the dataset."
    )
    parser.add_argument(
            "--seed",
            type=int,
            default=42
    )

    args = parser.parse_args()

    # test the model
    if args.data and args.num_vars:

        all_data : list[list[float]] = []

        # read data from file
        with open(args.data, "r") as file:
            for line in file:
                # ignore headers
                if line.strip().startswith("##"):
                    continue

                try:
                    row = [float(token) for token in line.split()]
                except ValueError as err:
                    sys.exit("*** Error: Data file not is expected format.\n")

                all_data.append(row)


        if not all_data:
            sys.exit(f"*** Error: No data found in {args.data}\n")


        all_data = np.asarray(all_data, dtype=float)
        var = all_data[:, :num_vars]
        resp = all_data[:, num_vars]

        # train-test split
        train_var, test_var, train_resp, test_resp = train_test_split(
                var,
                resp,
                test_size=0.2,
                random_state=args.seed,
                shuffle=True
        )

        # Create the GP surrogate
        gp_surr = GPSurrogate()

        # Construct
        gp_surr.construct(train_var, train_resp)

        # Loss
        gp_surr.loss(test_var, test_resp)

        print(f"Loss: {loss:16.8e}\n")

        # Save
        gp_surr.save() 


