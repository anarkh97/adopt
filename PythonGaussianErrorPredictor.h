/*  _______________________________________________________________________

    Dakota: Explore and predict with confidence.
    Copyright 2014-2024
    National Technology & Engineering Solutions of Sandia, LLC (NTESS).
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Dakota directory.
    _______________________________________________________________________ */

/************************************************************************
 * MODIFICATIONS IN THIS DERIVED FILE
 *
 * Derived from:
 *   - Dakota's dakota::surrogates::Python
 *   - Upstream commit: 0daaaa2237bd79445349c47d157a4a2b73db7452
 *
 * Modifications by Aditya Narkhede (AN), 2025:
 *   - Changed name from Python to PythonGaussianErrorPredictor
 *   - Removed dakota::surrogates namespace
 *   - Added a "variance" function to call the corresponding
 *     "predict_variance" from the python module.
 *   - Added "save_model" function to call the corresponding 
 *     "save" from python module. This is different from SurrogaBase 
 *     save as that function saves the PythonGaussianErrorPredictor
 *     object and not the underlying surrogate model.
 *
 * Additional Notes:
 *   - Keeping the original comments.
 *
 * Licensing:
 *   This file continues to be licensed under the GNU Lesser General Public
 *   License v2.1 or (at your option) any later version. See the bundled
 *   license file `LICENSE.LGPL-2.1` or
 *   <https://www.gnu.org/licenses/old-licenses/lgpl-2.1.html>.
 ************************************************************************/

#ifndef _PYTHON_GAUSSIAN_ERROR_PREDICTOR_H_
#define _PYTHON_GAUSSIAN_ERROR_PREDICTOR_H_

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/embed.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
namespace py = pybind11;

#include <SurrogatesBase.hpp>
#include <UtilDataScaler.hpp>
#include <UtilLinearSolvers.hpp>
#include <util_data_types.hpp>

/**
 *  \brief The PythonGaussianErrorPredictor class constructs a surrogate via python and has
 *         it ready for ADOPT  use.
 */

class PythonGaussianErrorPredictor : public dakota::surrogates::Surrogate {

 public:

  /**
   * \brief Constructor that sets moduleAndClassName and does not build.
   *
   * \param[in] module_and_class_name Name of python module file containing callback functions
   */
  PythonGaussianErrorPredictor(const std::string& module_and_class_name);

  /**
   * \brief Constructor sets moduleAndClassName and builds the python surrogate.
   *
   * \param[in] samples Matrix of data for surrogate construction - (num_samples
   * by num_features) \param[in] response Vector of targets for surrogate
   * construction - (num_samples by num_qoi = 1; only 1 response is supported
   * currently). \param[in] module_and_class_name Name of python module file
   * containing callback functions
   */
  PythonGaussianErrorPredictor(const Eigen::MatrixXd& samples, const Eigen::MatrixXd& response,
                       const std::string& module_and_class_name);

  /// Default destructor
  ~PythonGaussianErrorPredictor() { }

  /// Construct and populate the defaultConfigOptions.
  void default_options() override {
    return initialize_python();
  }

  // Allow disabling for field surrogates for now - RWH
  virtual bool diagnostics_available() override;

  /**
   * \brief Build the python surrogate using specified build data.
   *
   * \param[in] samples Matrix of data for surrogate construction - (num_samples
   * by num_features) \param[in] response Vector of targets for surrogate
   * construction - (num_samples by num_qoi = 1; only 1 response is supported
   * currently).
   */
  void build(const Eigen::MatrixXd& samples, const Eigen::MatrixXd& response) override;

  /**
   *  \brief Evaluate the scalar python surrogate at a set of prediction points.
   * \param[in] eval_points Matrix of prediction points - (num_pts
   * by num_features). \returns Values
   * of the python surrogate at the prediction points - (num_pts)
   */
  Eigen::VectorXd value(const Eigen::MatrixXd& eval_points) override;

  /**
   *  \brief Evaluate the field python surrogate at a set of prediction points.
   * \param[in] eval_points Matrix of prediction points - (num_pts
   * by num_features). \returns Values of the python surrogate at the
   * prediction points - (num_pts)
   */
  Eigen::VectorXd values(const Eigen::MatrixXd& eval_points) override;

  /**
   *  \brief Evaluate the gradient of the python scalar surrogate at a set of
   * prediction points. \param[in] eval_points Coordinates of
   * the prediction points - (num_pts by num_features). 
   * \returns Matrix of gradient
   * vectors at the prediction points - (num_pts by num_features).
   */
  Eigen::MatrixXd gradient(const Eigen::MatrixXd& eval_points) override;

  /**
   *  \brief Evaluate the Hessian of the python scalar surrogate at a single point.
   *  \param[in] eval_point Coordinates of the prediction point - (1 by
   * num_features). \returns Hessian matrix at the prediction point -
   *  (num_features by num_features).
   */
  Eigen::MatrixXd hessian(const Eigen::MatrixXd& eval_point) override;

  /**
   * Implemented by AN.
   */
  Eigen::VectorXd variance(const Eigen::MatrixXd& eval_points);

  double loss(const Eigen::MatrixXd& truth, const Eigen::MatrixXd& pred);

  void save_model(const std::string& filename);


  std::shared_ptr<dakota::surrogates::Surrogate> clone() const override {
    return std::make_shared<PythonGaussianErrorPredictor>(moduleAndClassName);
  }

 private:

  // --------------- PythonGaussianErrorPredictor Setup --------------------

  /// Name of python callback module file
  std::string moduleAndClassName;

  /// true if this class created the interpreter instance
  bool ownPython;

  /// true if python callback module is valid
  bool pyModuleActive;

  /// python Surrogate class
  py::object pySurrogate;

  /// flag for field-based surrogates
  bool isField;


  // -------------------------------------------------

  /// Initialize python interpreter and callback module
  void initialize_python();

  /// Verbosity level.
  int verbosity;
};

#endif
