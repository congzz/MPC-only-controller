windows_length = 20

class StateEstimator:
    def __init__(self, robot):
        self._robot = robot
        self._windows_length = windows_length
        self._linear_velocity_x_container = [0.0]*windows_length
        self._linear_velocity_y_container = [0.0]*windows_length
        self._linear_velocity_z_container = [0.0]*windows_length

        self._estimate_linear_velocity = [0.0]*3

    def run(self, current_time):
        linear_velocity_x, linear_velocity_y, linear_velocity_z = self._robot.getRobotLinearVelocity()

        del self._linear_velocity_x_container[0]
        del self._linear_velocity_y_container[0]
        del self._linear_velocity_z_container[0]

        self._linear_velocity_x_container.append(linear_velocity_x)
        self._linear_velocity_y_container.append(linear_velocity_y)
        self._linear_velocity_z_container.append(linear_velocity_z)

        self._estimate_linear_velocity[0] = sum(self._linear_velocity_x_container) / self._windows_length
        self._estimate_linear_velocity[1] = sum(self._linear_velocity_y_container) / self._windows_length
        self._estimate_linear_velocity[2] = sum(self._linear_velocity_z_container) / self._windows_length




