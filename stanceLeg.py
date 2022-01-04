import math
import sys

import numpy as np
import convexMPCController

import gait

mass = 108 / 9.8;
inertial = [0.01683993, 0.056579028, 0.064713601, 8.3902e-05, 0.000597679, 2.5134e-05]
planning_time_step = 0.025
weights = [5, 5, 0.2, 0.0, 0.0, 10, 0.0, 0.0, 1., 1., 1., 1., 0]
weight_alpha = 1e-5
friction_coeff = 0.4


class StanceLeg:
    def __init__(self, robot, gait_generator, state_estimator, ground_estimator):
        self._robot = robot
        self._ground_estimator = ground_estimator
        self._state_estimator = state_estimator
        self._gait_generator = gait_generator
        self._mpc_controller = convexMPCController.ContactForceMpcCaculate(mass, inertial,planning_time_step,
                                                        weights, weight_alpha, friction_coeff)

    def run(self, current_time, desired_x_speed, desired_y_speed, desired_twist_speed, desired_robot_height):
        leg_contact_state = [int(0)]*4
        for leg_id in range(4):
            if self._gait_generator._desired_leg_state[leg_id] == gait.STANCE:
                leg_contact_state[leg_id] = int(1)


        # print(leg_contact_state)
        # if leg_contact_state == [int(0)]*4:
        #     print(current_time)
        #     print(leg_contact_state, '\n\n')
        #     sys.exit(0)


        # linear_velocity_in_gravity_frame = self._robot.getRobotLinearVelocity()
        linear_velocity_in_gravity_frame = self._state_estimator._estimate_linear_velocity
        angular_velocity_in_gravity_frame = self._robot.getRobotAngulurVelocity()
        euler_angle = self._robot.getRobotEuler()
        euler_angle = list(euler_angle)
        euler_angle[2] = 0.0
        foot_pos_in_gravity_frame = self._robot.getFootPosInGravityFrame()
        foot_pos_in_gravity_frame = np.asarray(foot_pos_in_gravity_frame).flatten()
        robot_height = self._robot.getRobotPos()[2]

        desired_linear_velocity = self._ground_estimator._ground_posture_mat * np.mat([[desired_x_speed,
                                                                                        desired_y_speed,
                                                                                        0]]).transpose()
        desired_angular_velocity = [0.0, 0.0, desired_twist_speed]
        desired_euler_angle = [self._ground_estimator._ground_roll,
                               self._ground_estimator._ground_pitch,
                               0.0]

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
            if leg_contact_state[leg_id] == 0:
                continue

            self._robot.setFootForceInGravityFrame(leg_id, contact_force_trajectory[3*leg_id: 3*leg_id + 3])







