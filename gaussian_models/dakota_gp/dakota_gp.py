import sys
import argparse
from dakota import surrogates as daksurr
import numpy as np
from sklearn.model_selection import train_test_split
from typing import Dict, Any
from numpy.typing import ArrayLike 

class GPSurrogate:

    def __init__(self, params=None):
        nugget_options : Dict[str, Any] = {
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
        self.__gp = daksurr.GaussianProcess(self.config_options)

        return


    def construct(self, var : ArrayLike, resp : ArrayLike) -> None:

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


    def loss(self, truth : ArrayLike, pred : ArrayLike) -> float:

        loss : float = np.mean((truth - pred)**2)

        return np.sqrt(loss)

    
    def save(self, file_name : str) -> None:

        # Use Dakota's serialize output to save the GP model itself
        daksurr.save(self.__gp, file_name, False) # for binary output use True

        return


    def save_data(self, file_name : str, var : ArrayLike, resp : ArrayLike) -> None:

        num_points : int = var.shape[0]
        num_variables : int = var.shape[1]

        resp = resp.flatten()
        variance : ArrayLike = self.__gp.variance(var)
        pred_resp : ArrayLike = self.__gp.value(var)
        uncertainty : ArrayLike = np.sqrt(variance)
        loss : float = self.loss(var, resp)

        with open(file_name, "w") as file:
            # headers
            #file.write(f"## Iteration {self.__num_contruct_calls:04d}\n")
            file.write(f"## Loss {loss:16.8e}\n")
            file.write("## ")
            for i in range(num_variables):
                file.write(f"Parameter {i+1:04d}  |  ")
            file.write("NRMSE  |  Prediction  | Uncertainty\n")

            # data
            for i in range(num_points):
                for j in range(num_variables):
                    file.write(f"{var[i,j]:16.8e}  ")
                file.write(f"{resp[i]:16.8e}  {pred_resp[i]:16.8e}  {uncertainty[i]:16.8e}\n")

        return


# Simple code for testing ...
if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument(
            "--data",
            type=str,
            help="Data file containing variable and responses."
    )
    parser.add_argument(
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
        var = all_data[:, :args.num_vars]
        resp = all_data[:, args.num_vars]
        resp = resp.reshape(-1, 1) # responses should also be a 2D array

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

        # Predict
        test_pred = gp_surr.predict(test_var)

        # Loss
        loss : float = gp_surr.loss(test_resp, test_pred)

        print(f"Loss: {loss:16.8e}")

        # Save
        gp_surr.save("GaussianModel") 
        gp_surr.save_data("Train_data.txt", train_var, train_resp)
        gp_surr.save_data("Test_data.txt", test_var, test_resp)


