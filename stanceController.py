import numpy as np
import convexMPCController

mass = 108 / 9.8;
inertial = [0.01683993, 0.056579028, 0.064713601, 8.3902e-05, 0.000597679, 2.5134e-05]
planning_time_step = 0.03
weights = [5, 5, 0.2, 0.0, 0.0, 10, 0.0, 0.0, 1., 1., 1., 0., 0]
weight_alpha = 1e-5
friction_coeff = 0.4

class StanceController:
    def __init__(self, robot):
        self._robot = robot
        self._mpc_controller = convexMPCController.ContactForceMpcCaculate(mass, inertial, planning_time_step,
                                                                           weights, weight_alpha, friction_coeff)

    def run(self, desired_euler_angle, desired_robot_height, pos_x, pos_y, yaw):
        leg_contact_state = [int(1), int(1), int(1), int(1)]
        linear_velocity_in_gravity_frame = self._robot.getRobotLinearVelocity()
        angular_velocity_in_gravity_frame = self._robot.getRobotAngulurVelocity()
        euler_angle = self._robot.getRobotEuler()
        euler_angle = list(euler_angle)
        euler_angle[2] = 0.0
        foot_pos_in_gravity_frame = self._robot.getFootPosInGravityFrame()
        foot_pos_in_gravity_frame = np.asarray(foot_pos_in_gravity_frame).flatten()
        robot_height = self._robot.getRobotPos()[2]

        kp = 10
        desired_angular_velocity = [0.0]*3
        desired_angular_velocity[2] = kp * (yaw - self._robot.getRobotEuler()[2])

        desired_linear_velocity = [0.0]*3
        desired_linear_velocity[0] = kp * (pos_x - linear_velocity_in_gravity_frame[0])
        desired_linear_velocity[1] = kp * (pos_y - linear_velocity_in_gravity_frame[1])


        contact_force_trajectory = self._mpc_controller.getContactForceTrajectory(linear_velocity_in_gravity_frame,
                                                                                  angular_velocity_in_gravity_frame,
                                                                                  euler_angle,
                                                                                  foot_pos_in_gravity_frame,
                                                                                  leg_contact_state,
                                                                                  desired_linear_velocity,
                                                                                  desired_angular_velocity,
                                                                                  desired_euler_angle,
                                                                                  desired_robot_height)

        for leg_id in range(4):
            self._robot.setFootForceInGravityFrame(leg_id, contact_force_trajectory[3 * leg_id: 3 * leg_id + 3])
