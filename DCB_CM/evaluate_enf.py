

from common_evaluation import *
common_functions = Common_functions()


class Evaluate_ENF():
    def __init__(self, force: [], displacement: [], time_utm: [], time_delam: [], delam_length: [], spec_data: 'class', names_par: dict, name: str ):
        self.force = force
        self.displacement = displacement
        self.time_utm = time_utm
        self.time_delam = time_delam
        self.delam_length = delam_length
        self.spec_data = spec_data
        self.names_par = names_par
        self.name = name
        # add initial delamination to all values in delam_length vector
        self.delamination_expand()

    # ~~~~~~~~~~~~~~~~~~~~~~~~~
    # possible common functions
    #
    def Return_dimensions(self, name_dim: str) -> []:
    # Function returns list of dimensions specified by name in dictionary 
        name_index = self.names_par[name_dim]

        num_index = [self.spec_data.data_name.index(i) for i in name_index ]
        values = [self.spec_data.data_value[i]*self.spec_data.data_dimension[i] for i in num_index ]

        return values
    #
    def delamination_expand(self):
        Delam_values = self.Return_dimensions('delam')
        delam_a0 = min(Delam_values)
    
        # add intitial delamination to values
        Del_length_compl = [float(delam_a0) + float(x) for x in self.delam_length]
        self.delam_length = Del_length_compl
    # ~~~~~~~~~~~~~~~~~~~~~~~~~

    def compute_compliance_ENF(self, path: str, num_iter: [int]) -> float:
        # Computing compliance m for ENF 
    
        a = []
        C = []
    
        for counter, j in enumerate(num_iter):
            a.append(self.delam_length[counter]**3)
            C.append(self.displacement[j] / self.force[j])
    
        # c0 + c1*x + c2*x**2
        c0, c1 = np.polynomial.polynomial.Polynomial.fit(C, a, 1, domain=[])
        yFitted = np.polyval([c1, c0], C)
    
        delta_y = abs(yFitted[1] - yFitted[0])
        delta_x = abs(C[1] - C[0])
        m = delta_x / delta_y
    
        delta = -c0/c1
    
        # expanding the polyfit curve to be visible in plot
        x_axis = np.linspace(min(C), max(C), 200)
        y_fitted_long = c1*x_axis + c0 
    
        fig = plt.figure()
        ax = fig.add_subplot(111)
        plt.plot(x_axis, y_fitted_long)
        plt.plot(C, a, 'o')
    
        plt.legend(["PolyFit", "Experiments"])
        plt.ylabel('$a^3$ [N/m]', fontsize=12)
        plt.xlabel('$C$ [m]', fontsize=12)
        plt.title('ENF - Compliance', fontsize=14)
        plt.grid(True)
        plt.tight_layout()

        path_fig_png = path + "_ENF" + ".png"
        path_fig_eps = path + "_ENF" + ".eps"
    
    #    fig.savefig(path_fig_eps)
        fig.savefig(path_fig_png)
    
        return m
    
    
    # ----------------------------------------
    def vypocet_energie_ENF(self, m: float):
    
        Width_values = self.Return_dimensions('width')
        Delam_values = self.Return_dimensions('delam')
        B = min(Width_values)
        a0 = min(Delam_values)

        F_max = max(self.force)
    
        defor_energ = (3*m*(F_max)**2 *(a0)**2)/(2*B) # Norm D7905
    
        # Bernardin Thesis approach Un-known source
        # index_F_max = (np.argmax(F))
        # float_lst = [float(x) for x in epsilon_n]
        # deltaC = float_lst[index_F_max]  # posunuti hodnoty v bode F_max v [m]
        # defor_energ = (9*a_del**2 * F_max * deltaC)/(2*B*(2*L**3 + 3*a_del**3)) # Bernardin Thesis
    
        return defor_energ
    
    # ----------------------------------------
    
    
    def vypocet_energie_ENF_R(self, F_max: float, m: float, del_length: float):

        Width_values = self.Return_dimensions('width')
        B = min(Width_values)
    
        defor_energ = (3*m*(F_max)**2 *(del_length)**2)/(2*B) # Norm D7905
    
        # Bernardin Thesis approach Un-known source
        # index_F_max = (np.argmax(F))
        # float_lst = [float(x) for x in epsilon_n]
        # deltaC = float_lst[index_F_max]  # posunuti hodnoty v bode F_max v [m]
        # defor_energ = (9*a_del**2 * F_max * deltaC)/(2*B*(2*L**3 + 3*a_del**3)) # Bernardin Thesis
    
        return defor_energ
    
    # ----------------------------------------
    
    def create_R_cuve_ENF(self, time_delay_R, path: str, voltage: list, voltage_list):
    
        print(f'Evaluating ENF Fracture energie')
    
    #    num_iter_r = found_indexes(Time_exp_ar, Time_R_ar)
        num_iter_r_v = common_functions.found_indexes_Voltage(voltage)
    
        # indexes are striped to evaluated values from photos
        num_iter_r_v = num_iter_r_v[voltage_list[0]:(voltage_list[0]+voltage_list[1])]
    
        ENF_G_II = []
        delam_set = []
    
        F_max = max(self.force)
    
        # compute the compliance for one sample
        m = self.compute_compliance_ENF(path, num_iter_r_v)
        print('Compliance m for ENF is: ', m)
        G_II = self.vypocet_energie_ENF(m)
        print(f"Computed mode II fracture toughness according to norm D7905 is: {G_II}")

    
        # loop thorough every delamination length value
        for counter, index in enumerate(num_iter_r_v):
    
            force_del = self.force[index]
            epsilon_del = self.displacement[index]
            Delam_length = self.delam_length[counter]
            
#            logger.info("index " + str(index) + " force " + str(force_del) + " epsilon " + str(epsilon_del) + " time " + str(self.time_utm[index]))
    
            #G_II = vypocet_energie_ENF(F_max, a, m)
            G_II = self.vypocet_energie_ENF_R(force_del, m, Delam_length)
    
            ENF_G_II.append(G_II)
    
    
        return ENF_G_II, self.delam_length 
    

# ----------------------------------------
# Graphical evaluation
def graph_ENF_R_curves(F: [[]], epsilon: [[]], name_TAR, color_set, marker_set, line_set, path_fig, range_min):

    fig = plt.figure()
    ax = fig.add_subplot(111)

    # constant average value from range_x1
    avg_x_const, avg_y_const = average_curve(epsilon, F, range_min)
    # average function
    avg_x, avg_y = average_curve(epsilon, F)

    # plot all R curves
    for cons in range(0, len(F)):
        plt.plot(epsilon[cons], F[cons], color=color_set[cons], marker=marker_set[cons], linestyle=line_set[cons])

    plt.plot(avg_x, avg_y, color='green')
    
    plt.xlabel(r'Delamination legth $a$ [m]')
    plt.ylabel(r'Fracture toughness $G_II$ [J$m^{-2}$]')

    plt.legend(name_TAR)

    path_fig_png = path_fig + ".png"
    path_fig_eps = path_fig + ".eps"

    # Plotting the constant curve which is the mean value from set range
    plt.plot(avg_x_const, avg_y_const, linewidth=3, color='k', linestyle='--')
    print("graph_ENF_R was created")

    fig.savefig(path_fig_eps)
    fig.savefig(path_fig_png)






