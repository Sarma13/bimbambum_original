""""
*  __       __          --
* |__)||\/||__) /\ |\/||__)/  \|\/|
* |__)||  ||__)/--\|  ||__)\__/|  |
*
* This file is part of BIMBAMBUM.
*
* classify_domain.py
*
* Description:
*
*
* -----------------------------------------------------------------------------
* Copyright (C) 2023 Armand Sieber
*
* This program is free software: you can redistribute it and/or modify it
* under the terms of the GNU General Public License as published by the Free
* Software Foundation, either version 3 of the License, or (at your option)
* any later version.
*
* This program is distributed in the hope that it will be useful, but WITHOUT
* ANY WARRANTY; without even the implied warranty of  MERCHANTABILITY or
* FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
* more details.
*
* You should have received a copy of the GNU General Public License along with
* this program.  If not, see <http://www.gnu.org/licenses/>.
* -----------------------------------------------------------------------------
*
"""

import flow_field
import inputs

if __name__ == '__main__':

    init_data = inputs.simulation_inputs
    t_steps = init_data.time_step

    for i in range(len(t_steps)):
        t_i = t_steps[i]
        print('Time step: ', t_i)
        simulation = flow_field.flow_field(t_i)
        simulation.domain_creation()
        simulation.compute_velocities()
        simulation.compute_pressure()
        simulation.write_solution()
        simulation.post_processing()
