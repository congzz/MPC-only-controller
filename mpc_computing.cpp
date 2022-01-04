#include <cmath>
#include <Eigen>
#include <vector>
#include <qpOASES.hpp>

#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

const int planning_horizion_ = 10;

typedef Eigen::Matrix<qpOASES::real_t, Eigen::Dynamic, Eigen::Dynamic,
                      Eigen::RowMajor>
    RowMajorMatrixXd;


Eigen::Matrix<double, 13, 13> calculateGMat(Eigen::Matrix<double, 13, 13> A_mat, double time_step)
{
    Eigen::Matrix<double, 13, 13> G_mat;
    G_mat.setZero();

    Eigen::Matrix<double, 13, 13> A_mat_t = A_mat * time_step;

    Eigen::Matrix<double, 13, 13> term = Eigen::MatrixXd::Identity(13, 13);

    for(int i = 0; i < 20; i++)
    {
        G_mat += term;
        term *= 1 / (i + 1) * A_mat_t;
    }

    return G_mat;
}

Eigen::Matrix<double, 13, 12> calculateHMat(Eigen::Matrix<double, 13, 13> A_mat, 
                Eigen::Matrix<double, 13, 12> B_mat, double time_step)
{
        Eigen::Matrix<double, 13, 13> A_mat_t_int;
        A_mat_t_int.setZero();

        Eigen::Matrix<double, 13, 13> A_mat_t = A_mat * time_step;
        Eigen::Matrix<double, 13, 13> term = Eigen::MatrixXd::Identity(13, 13) * time_step;

        for(int i = 0; i < 20; i++)
        {
            A_mat_t_int += term;
            term *= 1 / (i + 2) * A_mat_t;
        }

        return A_mat_t_int * B_mat;
}

// Eigen::MatrixXd powMat(Eigen::MatrixXd mat, int n)
// {
//     Eigen::MatrixXd pow_mat = Eigen::MatrixXd::Identity(mat.rows(), mat.cols());

//     for (int i = 0; i < n; i++)
//     {
//         pow_mat = pow_mat * mat;
//     }

//     return pow_mat;
// }

void copyToMatrix(const Eigen::MatrixXd &input,
                  const std::vector<int> foot_contact_states,
                  int row_blk_size, int col_blk_size,
                  bool is_block_diagonal, Eigen::Map<RowMajorMatrixXd> *out)
{
    // the block index in the destination matrix.
    int row_blk = 0;
    for (int i = 0; i < planning_horizion_ * 4; ++i)
    {
        int leg_id = (i % 4);
        if (foot_contact_states[leg_id] == 0)
        {
            // skip the row block.
            continue;
        }
        if (is_block_diagonal)
        {
            // just copy the block
            int col_blk = row_blk;
            out->block(row_blk * row_blk_size, col_blk * col_blk_size, row_blk_size,
                       col_blk_size) = input.block(i * row_blk_size, i * col_blk_size,
                                                   row_blk_size, col_blk_size);
        }
        else
        {
            int col_blk = 0;
            // Non-diagonal, need to copy all elements.
            for (int j = 0; j < planning_horizion_ * 4; ++j)
            {
                int leg_id = (j % 4);
                if (foot_contact_states[leg_id] == 0)
                {
                    // skip the col block.
                    continue;
                }
                out->block(row_blk * row_blk_size, col_blk * col_blk_size, row_blk_size,
                           col_blk_size) =
                    input.block(i * row_blk_size, j * col_blk_size, row_blk_size,
                                col_blk_size);
                ++col_blk;
            }
        }
        ++row_blk;
    }
}

void copyToVec(const Eigen::VectorXd &vec,
               const std::vector<int> foot_contact_states,
               int blk_size,
               std::vector<qpOASES::real_t> *out)
{
    int buffer_index = 0;
    for (int i = 0; i < 4 * planning_horizion_; ++i)
    {
        int leg_id = (i % 4);
        if (foot_contact_states[leg_id] == 0)
        {
            // skip the block.
            continue;
        }
        // otherwise copy this block.

        for (int j = 0; j < blk_size; ++j)
        {
            int index = i * blk_size + j;
            (*out)[buffer_index] = vec[index];
            ++buffer_index;
        }
    }
}

Eigen::Matrix3d inertialMatAssign(std::vector<double> inertial)
{
    Eigen::Matrix3d inertial_mat;

    inertial_mat(0, 0) = inertial[0];
    inertial_mat(1, 1) = inertial[1];
    inertial_mat(2, 2) = inertial[2];

    inertial_mat(0, 1) = inertial_mat(1, 0) = inertial[3];
    inertial_mat(0, 2) = inertial_mat(2, 0) = inertial[4];
    inertial_mat(1, 2) = inertial_mat(2, 1) = inertial[5];

    return inertial_mat;
}

Eigen::Matrix3d constructMatMapFromAngularVelocityToEulerDot(std::vector<double> euler_angle)
{
    double cos_y = cos(euler_angle[1]);
    double tan_y = tan(euler_angle[1]);
    double cos_z = cos(euler_angle[2]);
    double sin_z = sin(euler_angle[2]);

    Eigen::Matrix3d mat_map;

    mat_map << cos_z / cos_y, sin_z / cos_y, 0,
        -sin_z, cos_z, 0,
        cos_z * tan_y, sin_z * tan_y, 1;

    return mat_map;
}

Eigen::Matrix3d constuctSkewMatFromVec(Eigen::Vector3d vec)
{
    Eigen::Matrix3d skew_mat;

    skew_mat << 0.0, -vec(2), vec(1),
        vec(2), 0.0, -vec(0),
        -vec(1), vec(0), 0.0;

    return skew_mat;
}

Eigen::Matrix3d rotx(double x)
{
    Eigen::Matrix3d rot_mat;

    double c = cos(x);
    double s = sin(x);

    rot_mat << 1.0, 0.0, 0.0,
        0.0, c, -s,
        0.0, s, c;

    return rot_mat;
}

Eigen::Matrix3d roty(double y)
{
    Eigen::Matrix3d rot_mat;

    double c = cos(y);
    double s = sin(y);

    rot_mat << c, 0.0, s,
        0.0, 1.0, 0.0,
        -s, 0.0, c;

    return rot_mat;
}

Eigen::Matrix3d rotz(double z)
{
    Eigen::Matrix3d rot_mat;

    double c = cos(z);
    double s = sin(z);

    rot_mat << c, -s, 0.0,
        s, c, 0.0,
        0.0, 0.0, 1.0;

    return rot_mat;
}

Eigen::Matrix3d convertInertialMatFromBodyFrameToInvGravityFrame(Eigen::Matrix3d inertial_in_body_frame,
                                                                 std::vector<double> euler_angle)
{
    Eigen::Matrix3d inertial_in_gravity_frame;

    Eigen::Matrix3d rot_zyx = rotz(euler_angle[2]) * roty(euler_angle[1]) * rotx(euler_angle[0]);

    inertial_in_gravity_frame = rot_zyx * inertial_in_body_frame * rot_zyx.transpose();

    Eigen::Matrix3d inv_inertial_in_gravity_frame = inertial_in_gravity_frame.inverse();

    return inv_inertial_in_gravity_frame;
}

class ContactForceMpcCaculate
{
public:
    ContactForceMpcCaculate(double mass, std::vector<double> inertial_in_body_frame,
                            double planning_time_step,
                            std::vector<double> weights, double weight_alpha, double friction_coeff);

    void updateState(std::vector<double> &linear_velocity_in_gravity_frame,
                     std::vector<double> &angular_velocity_in_gravity_frame,
                     std::vector<double> &euler_angle,
                     std::vector<double> &foot_pos_in_gravity_frame,
                     std::vector<int> &leg_contact_state,
                     std::vector<double> &desired_linear_velocity_in_gravity_frame,
                     std::vector<double> &desired_angular_velocity_in_gravity_frame,
                     std::vector<double> &desired_euler_angle,
                     double desired_robot_height);

    std::vector<double> getContactForceTrajectory(std::vector<double> linear_velocity_in_gravity_frame,
                                                  std::vector<double> angular_velocity_in_gravity_frame,
                                                  std::vector<double> euler_angle,
                                                  std::vector<double> foot_pos_in_gravity_frame,
                                                  std::vector<int> leg_contact_state,
                                                  std::vector<double> desired_linear_velocity_in_gravity_frame,
                                                  std::vector<double> desired_angular_velocity_in_gravity_frame,
                                                  std::vector<double> desired_euler_angle,
                                                  double desired_robot_height);
    void constructABMat();
    void constructGHMat();
    void constructAqpBqp();
    void constructReferenceTrajectory();

    void constructInitialState();
    void constructQPFormation();
    void constructQuadTermQP();

    // void constructGMatSparse();
    // void constructHMatSparse();

    void updateConstraintsMatrix();
    void calculateConstraintBounds();

    void calculateRobotHeight();

private:
    double inv_mass_;
    double mass_;
    Eigen::Matrix3d inertial_in_body_frame_;
    double planning_time_step_;
    Eigen::Matrix<double, 13 * planning_horizion_, 13 * planning_horizion_> weights_mat_;
    Eigen::Matrix<double, 13, 13> weights_single_mat_;
    Eigen::Matrix<double, 12 * planning_horizion_, 12 * planning_horizion_> weight_alpha_mat_;
    Eigen::Matrix<double, 12, 12> weight_alpha_single_mat_;
    Eigen::Matrix<double, 13, 13> A_mat_;
    Eigen::Matrix<double, 13, 12> B_mat_;
    Eigen::Matrix<double, 13, 13> G_mat_;
    Eigen::Matrix<double, 13, 12> H_mat_;
    Eigen::Matrix<double, 13 * planning_horizion_, 13> Aqp_;
    Eigen::Matrix<double, 13 * planning_horizion_, 12 * planning_horizion_> Bqp_;
    Eigen::Matrix<double, 13 * planning_horizion_, 1> reference_trajectory_;
    Eigen::Matrix<double, 13, 1> initial_state_;
    Eigen::Matrix<double, 12 * planning_horizion_, 12 * planning_horizion_> quad_term_qp_;
    Eigen::Matrix<double, 12 * planning_horizion_, 1> linear_term_qp_;
    std::vector<double> linear_velocity_in_gravity_frame_;
    std::vector<double> angular_velocity_in_gravity_frame_;
    std::vector<double> euler_angle_;
    std::vector<double> foot_pos_in_gravity_frame_;
    std::vector<int> leg_contact_state_;
    double robot_height_;
    std::vector<double> desired_linear_velocity_in_gravity_frame_;
    std::vector<double> desired_angular_velocity_in_gravity_frame_;
    std::vector<double> desired_euler_angle_;
    double desired_robot_height_;

    Eigen::Matrix<double, 4 * 5 * planning_horizion_, 12 * planning_horizion_> constraint_;
    Eigen::Matrix<double, 4 * 5 * planning_horizion_, 1 >constraint_lb_;
    Eigen::Matrix<double, 4 * 5 * planning_horizion_, 1 >constraint_ub_;

    Eigen::SparseMatrix<double> G_mat_sparse_;
    Eigen::SparseMatrix<double> H_mat_sparse_;

    std::vector<double> contact_force_trajectory_;

    double friction_coeff_;
};

ContactForceMpcCaculate::ContactForceMpcCaculate(double mass, std::vector<double> inertial_in_body_frame,
                                                 double planning_time_step,
                                                 std::vector<double> weights, double weight_alpha, double friction_coeff)
    : inv_mass_(1.0 / mass),
      mass_(mass),
      planning_time_step_(planning_time_step),
      G_mat_sparse_(13, 13),
      H_mat_sparse_(13, 12),
      contact_force_trajectory_(12 * planning_horizion_, 0),
      friction_coeff_(friction_coeff),
      robot_height_(0.0)
{
    inertial_in_body_frame_ = inertialMatAssign(inertial_in_body_frame);

    A_mat_.setZero();
    B_mat_.setZero();
    G_mat_.setZero();
    H_mat_.setZero();
    Aqp_.setZero();
    Bqp_.setZero();
    reference_trajectory_.setZero();
    quad_term_qp_.setZero();
    linear_term_qp_.setZero();
    constraint_.setZero();
    constraint_lb_.setZero();
    constraint_ub_.setZero();

    weights_mat_.setZero();
    weights_single_mat_.setZero();
    weight_alpha_mat_.setZero();
    weight_alpha_single_mat_.setZero();

    for (int i = 0; i < 13; i++)
    {
        weights_single_mat_(i, i) = weights[i];
    }

    weight_alpha_single_mat_ = Eigen::MatrixXd::Identity(12, 12) * weight_alpha;

    for (int i = 0; i < planning_horizion_; i++)
    {
        weights_mat_.block(13 * i, 13 * i, 13, 13) = weights_single_mat_;
        weight_alpha_mat_.block(12 * i, 12 * i, 12, 12) = weight_alpha_single_mat_;
    }
}

void ContactForceMpcCaculate::updateConstraintsMatrix()
 {
     Eigen::Matrix<double, 5, 3> block_mat_single_constraint;

     block_mat_single_constraint << 1.0, 0.0, -friction_coeff_,
        1.0, 0.0, friction_coeff_,
        0.0, 1.0, -friction_coeff_,
        0.0, 1.0, friction_coeff_,
        0.0, 0.0, 1.0; 
    for(int i = 0; i < 4 * planning_horizion_; i++)
    {
        constraint_.block(5*i,3*i, 5, 3) = block_mat_single_constraint;

    }
}

void ContactForceMpcCaculate::calculateConstraintBounds()
{
    Eigen::Matrix<double, 5, 1> block_mat_single_constraint_lb;
    Eigen::Matrix<double, 5, 1> block_mat_single_constraint_ub;

    double lb = -10.0 * mass_ * 9.8;
    double ub = -lb;

    block_mat_single_constraint_lb << lb, 0.0, lb, 0.0, 0.0;
    block_mat_single_constraint_ub << 0.0, ub, 0.0, ub, ub;

    for(int i = 0; i < 4 * planning_horizion_; i++)
    {
        constraint_lb_.block(5*i, 0, 5, 1) = block_mat_single_constraint_lb;
        constraint_ub_.block(5*i, 0, 5, 1) = block_mat_single_constraint_ub;
    }

}

void ContactForceMpcCaculate::updateState(std::vector<double> &linear_velocity_in_gravity_frame,
                                          std::vector<double> &angular_velocity_in_gravity_frame,
                                          std::vector<double> &euler_angle,
                                          std::vector<double> &foot_pos_in_gravity_frame,
                                          std::vector<int> &leg_contact_state,
                                          std::vector<double> &desired_linear_velocity_in_gravity_frame,
                                          std::vector<double> &desired_angular_velocity_in_gravity_frame,
                                          std::vector<double> &desired_euler_angle,
                                          double desired_robot_height)
{
    linear_velocity_in_gravity_frame_ = linear_velocity_in_gravity_frame;
    angular_velocity_in_gravity_frame_ = angular_velocity_in_gravity_frame;
    euler_angle_ = euler_angle;
    foot_pos_in_gravity_frame_ = foot_pos_in_gravity_frame;
    leg_contact_state_ = leg_contact_state;
    desired_linear_velocity_in_gravity_frame_ = desired_linear_velocity_in_gravity_frame;
    desired_angular_velocity_in_gravity_frame_ = desired_angular_velocity_in_gravity_frame;
    desired_euler_angle_ = desired_euler_angle;
    desired_robot_height_ = desired_robot_height;
}

void ContactForceMpcCaculate::constructABMat()
{
    A_mat_.block(0, 6, 3, 3) = constructMatMapFromAngularVelocityToEulerDot(euler_angle_);

    A_mat_.block(3, 9, 3, 3) = Eigen::MatrixXd::Identity(3, 3);

    A_mat_(11, 12) = 1.0;

    Eigen::Vector3d foot_pos_in_gravity_frame_0, foot_pos_in_gravity_frame_1;
    Eigen::Vector3d foot_pos_in_gravity_frame_2, foot_pos_in_gravity_frame_3;

    for (int i = 0; i < 3; i++)
    {
        foot_pos_in_gravity_frame_0(i) = foot_pos_in_gravity_frame_[i];
        foot_pos_in_gravity_frame_1(i) = foot_pos_in_gravity_frame_[i + 3];
        foot_pos_in_gravity_frame_2(i) = foot_pos_in_gravity_frame_[i + 6];
        foot_pos_in_gravity_frame_3(i) = foot_pos_in_gravity_frame_[i + 9];
    }

    Eigen::Matrix3d foot_pos_in_gravity_frame_0_skew = constuctSkewMatFromVec(foot_pos_in_gravity_frame_0);
    Eigen::Matrix3d foot_pos_in_gravity_frame_1_skew = constuctSkewMatFromVec(foot_pos_in_gravity_frame_1);
    Eigen::Matrix3d foot_pos_in_gravity_frame_2_skew = constuctSkewMatFromVec(foot_pos_in_gravity_frame_2);
    Eigen::Matrix3d foot_pos_in_gravity_frame_3_skew = constuctSkewMatFromVec(foot_pos_in_gravity_frame_3);

    Eigen::Matrix3d inv_inertial_in_gravity_frame = convertInertialMatFromBodyFrameToInvGravityFrame(inertial_in_body_frame_,
                                                                                                     euler_angle_);

    B_mat_.block(6, 0, 3, 3) = inv_inertial_in_gravity_frame * foot_pos_in_gravity_frame_0_skew;
    B_mat_.block(6, 3, 3, 3) = inv_inertial_in_gravity_frame * foot_pos_in_gravity_frame_1_skew;
    B_mat_.block(6, 6, 3, 3) = inv_inertial_in_gravity_frame * foot_pos_in_gravity_frame_2_skew;
    B_mat_.block(6, 9, 3, 3) = inv_inertial_in_gravity_frame * foot_pos_in_gravity_frame_3_skew;

    for (int i = 0; i < 4; i++)
    {
        B_mat_.block(9, 3 * i, 3, 3) = Eigen::MatrixXd::Identity(3, 3) * inv_mass_;
    }
}

void ContactForceMpcCaculate::constructGHMat()
{
    G_mat_ = calculateGMat(A_mat_, planning_time_step_);
    H_mat_ = calculateHMat(A_mat_,  B_mat_, planning_time_step_);
    // G_mat_ = Eigen::MatrixXd::Identity(13, 13) + A_mat_ * planning_time_step_;
    // H_mat_ = B_mat_ * planning_time_step_;
}


// void ContactForceMpcCaculate::constructGMatSparse()
// {
//     std::vector<Eigen::Triplet<double>> triplets;

//     for (int i = 0; i < 13; i++)
//     {
//         triplets.emplace_back(i, i, G_mat_(i, i));
//     }

//     for (int row = 0; row < 3; row++)
//     {
//         for (int col = 6; col < 9; col++)
//         {
//             triplets.emplace_back(row, col, G_mat_(row, col));
//         }
//     }

//     for (int i = 0; i < 3; i++)
//     {
//         triplets.emplace_back(3 + i, 9 + i, G_mat_(3 + i, 9 + i));
//     }

//     triplets.emplace_back(11, 12, G_mat_(11, 12));

//     G_mat_sparse_.setFromTriplets(triplets.begin(), triplets.end());
// }

// void ContactForceMpcCaculate::constructHMatSparse()
// {
//     std::vector<Eigen::Triplet<double>> triplets;

//     for (int row = 6; row < 9; row++)
//     {
//         for (int col = 0; col < 12; col++)
//         {
//             triplets.emplace_back(row, col, H_mat_(row, col));
//         }
//     }

//     for (int i = 0; i < 4; i++)
//     {
//         for (int j = 0; j < 3; j++)
//         {
//             triplets.emplace_back(j + 9, 3 * i + j, H_mat_(j + 9, 3 * i + j));
//         }
//     }

//     H_mat_sparse_.setFromTriplets(triplets.begin(), triplets.end());
// }

void ContactForceMpcCaculate::constructAqpBqp()
{
    Eigen::Matrix<double, 13, 13> G_pow = G_mat_;


    for (int i = 0; i < planning_horizion_; i++)
    {
        Aqp_.block(i * 13, 0, 13, 13) = G_pow;
        G_pow = G_pow * G_mat_;
    }

 
    Eigen::Matrix<double, 13, 12> G_pow_H = H_mat_;

    for (int i = 0; i < planning_horizion_; i++)
    {
        for (int j = 0; j < planning_horizion_ - i; j++)
        {
            Bqp_.block(13 * i + 13 * j, 12 * j, 13, 12) = G_pow_H;
        }
        G_pow_H = G_mat_ * G_pow_H;
    }
}

void ContactForceMpcCaculate::constructReferenceTrajectory()
{
    Eigen::Vector3d desired_angular_velocity_in_gravity_frame_vec(desired_angular_velocity_in_gravity_frame_[0],
                                                                  desired_angular_velocity_in_gravity_frame_[1],
                                                                  desired_angular_velocity_in_gravity_frame_[2]);

    Eigen::Vector3d desired_euler_dot = constructMatMapFromAngularVelocityToEulerDot(desired_euler_angle_) *
                                        desired_angular_velocity_in_gravity_frame_vec;

    for (int i = 0; i < planning_horizion_; i++)
    {
        reference_trajectory_(13 * i) = desired_euler_angle_[0] + desired_euler_dot(0) * (i + 1) * planning_time_step_;
        reference_trajectory_(13 * i + 1) = desired_euler_angle_[1] + desired_euler_dot(1) * (i + 1) * planning_time_step_;
        reference_trajectory_(13 * i + 2) = desired_euler_angle_[2] + desired_euler_dot(2) * (i + 1) * planning_time_step_;

        reference_trajectory_(13 * i + 3) = desired_linear_velocity_in_gravity_frame_[0] * (i + 1) * planning_time_step_;
        reference_trajectory_(13 * i + 4) = desired_linear_velocity_in_gravity_frame_[1] * (i + 1) * planning_time_step_;
        reference_trajectory_(13 * i + 5) = desired_linear_velocity_in_gravity_frame_[2] * (i + 1) * planning_time_step_ + desired_robot_height_;

        reference_trajectory_(13 * i + 6) = desired_angular_velocity_in_gravity_frame_[0];
        reference_trajectory_(13 * i + 7) = desired_angular_velocity_in_gravity_frame_[1];
        reference_trajectory_(13 * i + 8) = desired_angular_velocity_in_gravity_frame_[2];

        reference_trajectory_(13 * i + 9) = desired_linear_velocity_in_gravity_frame_[0];
        reference_trajectory_(13 * i + 10) = desired_linear_velocity_in_gravity_frame_[1];
        reference_trajectory_(13 * i + 11) = desired_linear_velocity_in_gravity_frame_[2];

        reference_trajectory_(13 * i + 12) = -9.8;
    }
}

void ContactForceMpcCaculate::constructInitialState()
{
    initial_state_(0) = euler_angle_[0];
    initial_state_(1) = euler_angle_[1];
    initial_state_(2) = euler_angle_[2];

    initial_state_(3) = 0.0;
    initial_state_(4) = 0.0;
    initial_state_(5) = robot_height_;

    initial_state_(6) = angular_velocity_in_gravity_frame_[0];
    initial_state_(7) = angular_velocity_in_gravity_frame_[1];
    initial_state_(8) = angular_velocity_in_gravity_frame_[2];

    initial_state_(9) = linear_velocity_in_gravity_frame_[0];
    initial_state_(10) = linear_velocity_in_gravity_frame_[1];
    initial_state_(11) = linear_velocity_in_gravity_frame_[2];

    initial_state_(12) = -9.8;
}

void ContactForceMpcCaculate::constructQPFormation()
{
   

    Eigen::DiagonalMatrix<double, 13 * planning_horizion_> weights_diag_(weights_mat_.diagonal());
    linear_term_qp_ = Bqp_.transpose() * (weights_diag_ * (Aqp_ * initial_state_ - reference_trajectory_)); // 0.0001

    
    constructQuadTermQP(); // 0.0017
    
}

void ContactForceMpcCaculate::constructQuadTermQP()
{


    Eigen::DiagonalMatrix<double, 13> weights_single_diag_(weights_single_mat_.diagonal());
    Eigen::DiagonalMatrix<double, 12> weight_alpha_single_diag_(weight_alpha_single_mat_.diagonal());

    Eigen::Matrix<double, 13 * planning_horizion_, 12> G_pow_H = Bqp_.block(0, 0, 13 * planning_horizion_, 12);
    for (int i = planning_horizion_ - 1; i >= 0; --i)
    {
        quad_term_qp_.block(i * 12, (planning_horizion_ - 1) * 12, 12, 12) =
            G_pow_H.block((planning_horizion_ - i - 1) * 13, 0, 13, 12).transpose() *
            weights_single_diag_ * H_mat_;
        // Fill the lower-triangle part by transposing the corresponding
        // upper-triangle part.
        if (i != planning_horizion_ - 1)
        {
            quad_term_qp_.block((planning_horizion_ - 1) * 12, i * 12, 12,
                                12) =
                quad_term_qp_
                    .block(i * 12, (planning_horizion_ - 1) * 12, 12,
                           12)
                    .transpose();
        }
    }

    // We then fill in the submatrices in the middle by propagating the values
    // from lower right to upper left.
  
    for (int i = planning_horizion_ - 2; i >= 0; --i)
    {
        // Diagonal block.
        quad_term_qp_.block(i * 12, i * 12, 12, 12) =
            quad_term_qp_.block((i + 1) * 12, (i + 1) * 12, 12,
                                12) +
            G_pow_H.block((planning_horizion_ - i - 1) * 13, 0, 13, 12)
                    .transpose() *
                weights_single_diag_ *
                G_pow_H.block((planning_horizion_ - i - 1) * 13, 0, 13,
                              12);
        // Off diagonal block
        for (int j = i + 1; j < planning_horizion_ - 1; ++j)
        {
            quad_term_qp_.block(i * 12, j * 12, 12, 12) =
                quad_term_qp_.block((i + 1) * 12, (j + 1) * 12, 12,
                                    12) +
                G_pow_H.block((planning_horizion_ - i - 1) * 13, 0, 13, 12)
                        .transpose() *
                    weights_single_diag_ *
                    G_pow_H.block((planning_horizion_ - j - 1) * 13, 0, 13,
                                  12);
            // Fill the lower-triangle part by transposing the corresponding
            // upper-triangle part.
            quad_term_qp_.block(j * 12, i * 12, 12, 12) =
                quad_term_qp_.block(i * 12, j * 12, 12, 12)
                    .transpose();
        }
    }


    // Multiply by 2 and add alpha.
    for (int i = 0; i < planning_horizion_; ++i)
    {
        quad_term_qp_.block(i * 12, i * 12, 12, 12) +=
            weight_alpha_single_diag_;
    }

}

void ContactForceMpcCaculate::calculateRobotHeight()
{
    double sum = 0;
    int num = 0;
    for(int i = 0; i < 4; i++)
    {
        if(leg_contact_state_[i] == 1)
        {
            sum += foot_pos_in_gravity_frame_[3*i+2];
            num += 1;
        }
    }

    robot_height_ = fabs(sum / num);
}


std::vector<double> ContactForceMpcCaculate::getContactForceTrajectory(std::vector<double> linear_velocity_in_gravity_frame,
                                                                       std::vector<double> angular_velocity_in_gravity_frame,
                                                                       std::vector<double> euler_angle,
                                                                       std::vector<double> foot_pos_in_gravity_frame,
                                                                       std::vector<int> leg_contact_state,
                                                                       std::vector<double> desired_linear_velocity_in_gravity_frame,
                                                                       std::vector<double> desired_angular_velocity_in_gravity_frame,
                                                                       std::vector<double> desired_euler_angle,
                                                                       double desired_robot_height)
{

    updateState(linear_velocity_in_gravity_frame,
                angular_velocity_in_gravity_frame,
                euler_angle,
                foot_pos_in_gravity_frame,
                leg_contact_state,
                desired_linear_velocity_in_gravity_frame,
                desired_angular_velocity_in_gravity_frame,
                desired_euler_angle,
                desired_robot_height);

    constructABMat();
    constructGHMat();

    constructAqpBqp(); // 0.0001

    calculateRobotHeight();

    constructInitialState();
    constructReferenceTrajectory();

    constructQPFormation(); // 0.0003
    
    updateConstraintsMatrix();
    calculateConstraintBounds();

    

    int num_legs_in_contact = 0;
    for (auto &leg_contact : leg_contact_state)
    {
        num_legs_in_contact += leg_contact;
    }


    const int qp_dim = num_legs_in_contact * 3 * planning_horizion_;
    const int constraint_dim = num_legs_in_contact * 5 * planning_horizion_;

    std::vector<qpOASES::real_t> hessian(qp_dim * qp_dim, 0);
    Eigen::Map<RowMajorMatrixXd> hessian_mat_view(hessian.data(), qp_dim, qp_dim);
    // // Copy to the hessian
    copyToMatrix(quad_term_qp_, leg_contact_state, 3, 3, false, &hessian_mat_view);

    std::vector<qpOASES::real_t> g_vec(qp_dim, 0);
    // // Copy the g_vec
    copyToVec(linear_term_qp_, leg_contact_state, 3, &g_vec);

    std::vector<qpOASES::real_t> a_mat(qp_dim * constraint_dim, 0);
    Eigen::Map<RowMajorMatrixXd> a_mat_view(a_mat.data(), constraint_dim, qp_dim);
    copyToMatrix(constraint_, leg_contact_state, 5, 3, true, &a_mat_view);

    std::vector<qpOASES::real_t> a_lb(constraint_dim, 0);
    copyToVec(constraint_lb_, leg_contact_state, 5, &a_lb);

    std::vector<qpOASES::real_t> a_ub(constraint_dim, 0);
    copyToVec(constraint_ub_, leg_contact_state, 5, &a_ub);

    auto qp_problem = qpOASES::QProblem(qp_dim, constraint_dim, qpOASES::HST_UNKNOWN,
                               qpOASES::BT_TRUE);

    qpOASES::Options options;
    options.setToMPC();
    options.printLevel = qpOASES::PL_NONE;
    qp_problem.setOptions(options);

    
    int max_solver_iter = 100;

    // clock_t begin, end;

    // begin = clock();

    qp_problem.init(hessian.data(), g_vec.data(), a_mat.data(), nullptr,
                    nullptr, a_lb.data(), a_ub.data(), max_solver_iter,
                    nullptr);

    //  end = clock();

    // std::cout << double(end - begin) / CLOCKS_PER_SEC * 1000.0 << "\n";

    std::vector<qpOASES::real_t> qp_sol(qp_dim, 0);
    qp_problem.getPrimalSolution(qp_sol.data());
    for (auto &force : qp_sol)
    {
        force = -force;
    }

    int buffer_index = 0;
    for (int i = 0; i < 4 * planning_horizion_; ++i)
    {
        int leg_id = i % 4;
        if (leg_contact_state[leg_id] == 0)
        {
            contact_force_trajectory_[i * 3] = 0;
            contact_force_trajectory_[i * 3 + 1] = 0;
            contact_force_trajectory_[i * 3 + 2] = 0;
        }
        else
        {
            contact_force_trajectory_[i * 3] = qp_sol[buffer_index * 3];
            contact_force_trajectory_[i * 3 + 1] = qp_sol[buffer_index * 3 + 1];
            contact_force_trajectory_[i * 3 + 2] = qp_sol[buffer_index * 3 + 2];
            ++buffer_index;
        }
    }

    return contact_force_trajectory_;
}

PYBIND11_MODULE(convexMPCController, m) {
      
      
  py::class_<ContactForceMpcCaculate>(m, "ContactForceMpcCaculate")
      .def(py::init<double, std::vector<double>, double,
          std::vector<double>, double,double>())
      .def("getContactForceTrajectory", &ContactForceMpcCaculate::getContactForceTrajectory);
}
