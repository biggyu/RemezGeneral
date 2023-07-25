#include "RemezGeneral.h"
#include <NTL/RR.h>

RemezGeneral::RemezGeneral(RemezParam _params, long _section_num, double *_sections, long _deg)
    : params(_params), section_num(_section_num), sections(_sections), deg(_deg) {

    // for(int i = 0; i < _section_num * 2; i++) {
    //     cout << sections[i] << " ";
    // }
    // cout << endl;

    width = new RR[section_num];
    sc = new RR[section_num];

    //chebeval gives accurate result when boundary is [-1, 1].
    //as the boundary is different due to generalization, chebeval_k is defined to change the original range to the above.
    chebeval_k = 0;
    for (long i = 0; i < section_num; i++) {
        width[i] = abs(RR(sections[2 * i + 1] - sections[2 * i]));
        sc[i] = width[i] / pow(2.0, params.log_scan_step_diff);
        chebeval_k = abs(sections[2 * i]) > chebeval_k ? abs(sections[2 * i]) : chebeval_k;
        chebeval_k = abs(sections[2 * i + 1]) > chebeval_k ? abs(sections[2 * i + 1]) : chebeval_k;
    }
    chebeval_k = ceil(chebeval_k);
    cout << "chebeval_k: " << chebeval_k << endl;

    // for(int i = 0; i < section_num; i++) {
    //     cout << width[i] << " ";
    // }
    // cout << endl;

    // for(int i = 0; i < section_num; i++) {
    //     cout << sc[i] << " ";
    // }
    // cout << endl;

    approx = power2_RR(-params.log_approx_degree);

    max_err = 1000;
    min_err = 1;

    sample_point = new Point[deg + 2];
    extreme_point = new Point[2 * deg + section_num + 1];

    coeff = new RR[deg + 1];
}

void RemezGeneral::better_initialize() {
    int nodecount[section_num];
    int deg_bdd = deg + 2;

    for (long i = 0; i < section_num; i++) {
        nodecount[i] = 1;
    }
    int tot_deg = section_num;

    double *err = new double[section_num];
    for (long i = 0; i < section_num; i++) {
        err[i] = sections[2 * i + 1] - sections[2 * i];
    }

    RR::SetPrecision(params.RR_prec);

    double bdd[section_num];

    for (long i = 0; i < section_num; i++) {
        double temp = 0;
        for (long j = 1; j <= section_num; j++) {
            temp -= log2((double)j);
        }

        temp += section_num * log2(2 * M_PI);
        temp += log2(err[i]);

        bdd[i] = temp;
        for (int j = 1; j <= section_num - 1 - i; j++) {
            bdd[i] += log2((double)j + err[i]);
        }
        for (int j = 1; j <= section_num - 1 + i; j++) {
            bdd[i] += log2((double)j + err[i]);
        }
    }

    int max_iter = 200;
    // int iter;

    for (int iter = 0; iter < max_iter; iter++) {
        if (tot_deg >= deg_bdd) {
            break;
        }
        int maxi = max_index(bdd, section_num);

        if(maxi != 0) {
			if((tot_deg+2) > deg_bdd) {
				break;
            }
	
			for(int i = 0; i < section_num; i++) {
				bdd[i] -= log2(tot_deg + 1);
				bdd[i] -= log2(tot_deg + 2);
				bdd[i] += 2.0 * log2(2.0 * M_PI);

                bdd[i] += (i != maxi) ? (log2(abs((double)(i - maxi)) + err[i] * ((double)(i + maxi) + err[i]))) : 
                (log2(err[i] * 2.0 * (double)i + err[i]) - 1.0);
				// if(i != maxi) {	
				// 	bdd[i] += log2(abs((double)(i - maxi)) + err[i]);
				// 	bdd[i] += log2((double)(i + maxi) + err[i]);
				// } else { 
				// 	bdd[i] += (log2(err[i]) - 1.0);
				// 	bdd[i] += log2(2.0 * (double)i + err[i]);
				// }
			}
			
			tot_deg += 2;
		} else { 
			bdd[0] -= log2(tot_deg + 1);
			bdd[0] += (log2(err[0]) - 1.0);
			bdd[0] += log2(2.0 * M_PI);
			for(int i = 1; i < section_num; i++) {	
				bdd[i] -= log2(tot_deg + 1);
				bdd[i] += log2(2.0 * M_PI);
				bdd[i] += log2((double)i + err[i]);	
			}
		
			tot_deg += 1;	
		}
		
		nodecount[maxi] += 1;
    }
    // delete[] bdd;

    // cout << "node count" << endl;
    // for(int i = 0; i < (section_num + 1) / 2; i++) {
    //     cout << nodecount[i] << " ";
    // }
    // cout << endl;

    // cout << "sample pnt" << endl;
    // for(int i = 0; i < deg + 2; i++) {
    //     cout << sample_point[i].x << " " << sample_point[i].y << endl;
    // }
    // cout << endl;

    if (tot_deg == deg_bdd - 1) {
        nodecount[0]++;
        tot_deg++;
    }

    RR inter_size[section_num];
    for (long i = 0; i < section_num; i++) {
        inter_size[i] = width[i];
    }

    int cnt = 0;
    if ((nodecount[0] % 2) != 0) {
        sample_point[cnt++].x = RR(0);
    }

    for(int i = 0; i < section_num; i++) {
		for(int j = 1; j <= nodecount[i]; j++) {
			RR temp = ((RR(2 * j - 1)) * ComputePi_RR()) / (RR(2 * nodecount[i]));		
			sample_point[cnt++].x = RR(sections[2 * i]) + inter_size[i] * cos(temp);
			// sample_point[cnt++].x = RR(-i) - inter_size[i] * cos(temp);
		}
	}

    // for (int j = 1; j <= (nodecount[0] / 2); j++) {
    //     RR temp = ((RR(2 * j - 1)) * ComputePi_RR()) / (RR(2 * nodecount[0]));
    //     sample_point[cnt++].x = inter_size[i] * cos(temp);
    //     sample_point[cnt++].x = -inter_size[i] * cos(temp);
    //     // i = (i >= section_num) ? 0 : i + 1;
    // }

    sort(sample_point, sample_point + (deg + 2), xcompare);

    for (int i = 0; i < deg + 2; i++) {
        sample_point[i].y = function_value(sample_point[i].x);
        std::cout << sample_point[i].x  << std::endl;
    }

    // delete[] nodecount;
}

void RemezGeneral::initialize() {
	RR::SetPrecision(params.RR_prec);
	long* nodecount = new long[section_num];
	long ind = section_num - 1;
	bool alter = true;
	for(long j = 0; j < section_num; j++) {
		nodecount[j] = 0;
	}
	for(long j = 0; j < deg + 2; j++) {
		if(ind > 0) {
			nodecount[ind]++;
			ind--;
		} else {
			if(alter) {
				nodecount[ind]++;
				ind = section_num - 1;
				alter = false;
			}
			else {
				nodecount[section_num - 1]++;
				ind = section_num - 2;
				alter = true;
			}
		}
	}
	for(long j = 0; j < section_num; j++) {
		cout << nodecount[j] << " ";
	}
	cout << endl;
	ind = 0;
	for(int j = 0; j < section_num; j++) {
		for(long k = 0; k < nodecount[j]; k++) {
            sample_point[ind].x = sections[2 * j] + width[j] / (nodecount[j] + 1) * (k + 1);
			sample_point[ind].y = function_value(sample_point[ind].x);
			ind++;
		}
//		cout << endl;
	}

    cout << "sample points " << " ";
	for(long j = 0; j < deg + 2; j++) {
	   cout << sample_point[j].x << " ";
	}
    cout << endl;
}

void RemezGeneral::getcoeffwitherr() {
    RR::SetPrecision(params.RR_prec);

    vec_RR v, w, v_0;
    mat_RR m;

    v.SetLength(deg + 2);
    w.SetLength(deg + 2);
    v_0.SetLength(deg + 2);
    m.SetDims(deg + 2, deg + 2);

    for(long i = 0; i < deg + 2; i++) {
        v[i] = sample_point[i].x;
        w[i] = sample_point[i].y;
    }

    RR var;
    for(long i = 0; i < deg + 2; i++) {
        var = v[i];
        m[i][0] = 1;
        m[i][1] = var / chebeval_k;
        for (long j = 2; j < deg + 1; j++) {
            m[i][j] = 2 * m[i][1] * m[i][j - 1] - m[i][j - 2];
        }
        m[i][deg + 1] = 2 * (i % 2) - 1;
    }
    RR deter;
    mat_RR mtrns;
    mat_RR minv;
    transpose(mtrns, m);
    inv(deter, minv, mtrns);
    v_0 = w * minv;
    for (int i = 0; i <= deg; i++) {
        coeff[i] = v_0[i];
    }
    current_err = floor(power2_RR(params.log_round_prec) * abs(v_0[deg + 1])) / power2_RR(params.log_round_prec);
    cout << "coeffwitherr: " << current_err << endl;
}

void RemezGeneral::getextreme_local(Point *local_extreme_point, long &local_extreme_count, long section_ind) {
    RR::SetPrecision(params.RR_prec);
	RR::SetOutputPrecision(20);
	// cout << "getextreme start" << endl;
	long inc_1 = 0, inc_2 = 0;
	long tmpinc;
	RR scan_1, scan_2;
	
	scan_2 = sections[2 * section_ind];
	

	RR scan_y1, scan_y2 = chebeval(deg, coeff, scan_2 / chebeval_k) - function_value(scan_2);
	local_extreme_count = 0;
	
	RR detail[3];
	RR prec_sc, prec_ext, prec_x, tmp;
	long prec_iter, prec_ind;
	bool prec_end;
	long tmp_inc;

	while(scan_2 < sections[2 * section_ind + 1] + sc[section_ind]) {
		scan_1 = scan_2;
		scan_2 = scan_1 + sc[section_ind]; //move forward sc(epsilon)
		// std::cout << scan_2 << std::endl;


        inc_1 = inc_2;
		scan_y1 = scan_y2;
		scan_y2 = chebeval(deg, coeff, scan_2 / chebeval_k) - function_value(scan_2);
        
		// if(scan_y1 < scan_y2) {
        //     inc_2 = 1;
        // } else if(scan_y1 > scan_y2) {
        //     inc_2 = -1;
        // } else {
        //     inc_2 = 0;
        // }
        inc_2 = (scan_y1 < scan_y2) ? 1 : ((scan_y1 > scan_y2) ? -1 : 0);

		// cout << "inc : " << inc_2 << endl;
		// if((inc_1 == 1 && inc_2 != 1) || (inc_1 == -1 && inc_2 != -1) || inc_1 == 0) {
        if(inc_1 * inc_2 != 1) {
			prec_end = false;
			tmp_inc = inc_2;
			prec_x = scan_2;
			while(!prec_end) {
				prec_sc = (prec_x - scan_1) / 2;
				prec_end = true;
				if(inc_1 != 0) {
					// cout << "extreme detected" << endl;
					for(long k = 0; k < 3; k++) {
						detail[k] = scan_1 + (k - 1) * prec_sc;
					}
				} else {
					// cout << "boundary start" << endl;
					for(long k = 0; k < 3; k++) {
						detail[k] = scan_1 + prec_sc * k;
					}
				}
				prec_iter = 0;
				while(prec_iter < params.binary_prec) {
					prec_ext = chebeval(deg, coeff, detail[0] / chebeval_k) - function_value(detail[0]);
					prec_ind = 0;
					for(long k = 1; k < 3; k++) {
						tmp = chebeval(deg, coeff, detail[k] / chebeval_k) - function_value(detail[k]);
						if((inc_1 == 1 && prec_ext < tmp) || (inc_1 == -1 && prec_ext > tmp)) {
							prec_ext = tmp;
							prec_ind = k;
						} else if(inc_1 == 0) {
							if((inc_2 == 1 && prec_ext > tmp) || (inc_2 == -1 && prec_ext < tmp)) {
								prec_ext = tmp;
								prec_ind = k;
							}
						}
					}
					if(inc_1 == 0 && prec_ind != 0) {
                        prec_end = false;
                    }
					prec_x = detail[prec_ind];
					prec_sc = prec_sc / 2;
					// cout << "prec_ind : " << prec_ind << " prec_x : " << prec_x << endl;
					for(long k = 0; k < 3; k++) {
						// if(inc_1 == 0) {
                        //     detail[k] = max(prec_x - prec_sc, scan_1) + prec_sc * k;
                        // } else {
                        //     detail[k] = prec_x - prec_sc + prec_sc * k;
                        // }
                        detail[k] = (inc_1 == 0) ? max(prec_x - prec_sc, scan_1) + prec_sc * k : prec_x + (k - 1) * prec_sc;
					}
					prec_iter++;
				}
				// if(inc_2 == 1) {
                //     tmpinc = -1;
                // } else {
                //     tmpinc = 1;
                // }
                tmpinc = (inc_2 == 1) ? -1 : 1;
                
				// cout << "###" << inc_2 << " " << prec_ext << " " << tmpinc * prec_ext << endl;
				if(tmpinc * prec_ext >= current_err) {
					local_extreme_point[local_extreme_count].x = prec_x;
					local_extreme_point[local_extreme_count].y = prec_ext;
					// if(inc_2 == 1) {
                    //     local_extreme_point[local_extreme_count].locmm = -1;
                    // } else {
                    //     local_extreme_point[local_extreme_count].locmm = +1;
                    // }
                    local_extreme_point[local_extreme_count++].locmm = (inc_2 == 1) ? -1 : 1;

					// cout << extreme_point[local_extreme_count].x << " " << extreme_point[local_extreme_count].y << " " << extreme_point[local_extreme_count].locmm << endl;
					// local_extreme_count++;
					// cout << "prec_end : " << prec_end << endl;

                    cout << "[" << section_ind << "]: " << local_extreme_count << endl;
				}
				// if(!prec_end) { 
				// 	inc_2 = (-1) * inc_2; 
				// 	// cout << "inc : " << inc_2 << endl;
				// }
                inc_2 *= prec_end ? 1 : -1;
			}
			inc_2 = tmp_inc;
		}

        //밑에 조건문은 극값이 이전 step과 경계값 사이에 있을 경우에 극값을 찾기 위해 구현한 코드임. 이는 더 큰 sc, 즉 더 구간을 더 자세히 찾는 것으로 대체하여 주석처리함.
		// if(fracpart(scan_2) > width[section_ind] + sc[section_ind] / 2) { ///////////
		// 	// cout << "boundary detected" << endl;
		// 	// cout << "inc : " << inc_2 << endl;
		// 	scan_1 = round(scan_1) + width[section_ind];
		// 	scan_2 = round(scan_2) + 1 - width[section_ind];
		// 	scan_y1 = chebeval(deg, coeff, scan_1 / chebeval_k) - function_value(scan_1); ///////
		// 	scan_y2 = chebeval(deg, coeff, scan_2 / chebeval_k) - function_value(scan_2); ///////

		// 	prec_end = false;
		// 	prec_x = scan_1 - sc[section_ind];
		// 	while(!prec_end) {
		// 		prec_sc = (scan_1 - prec_x) / 2;
		// 		prec_end = true;
		// 		for(long k = 0; k < 3; k++) {
		// 			detail[k] = scan_1 + (k - 2) * prec_sc;
		// 		}
		// 		prec_iter = 0;
		// 		while(prec_iter < params.binary_prec) { // < 10
		// 			prec_ext = chebeval(deg, coeff, detail[0] / chebeval_k) - function_value(detail[0]);
		// 			prec_ind = 0;
		// 			for(long k = 1; k < 3; k++) {
		// 				tmp = chebeval(deg, coeff, detail[k] / chebeval_k) - function_value(detail[k]);
		// 				if((inc_2 == 1 && prec_ext < tmp) || (inc_2 == -1 && prec_ext > tmp)) {
		// 					prec_ext = tmp;
		// 					prec_ind = k;
		// 				}
		// 			}
		// 			if(prec_ind != 2) {
        //                 prec_end = false;
        //             }
		// 			prec_x = detail[prec_ind];
		// 			prec_sc /= 2;
		// 			for(long k = 0; k < 3; k++) {
		// 				detail[k] = min(prec_x + prec_sc, scan_1) + (k - 2) * prec_sc;
		// 			}
		// 			prec_iter++;
		// 			// cout << "prec_ind : " << prec_ind << " prec_x : " << prec_x << endl;
		// 		}
		// 		// if(inc_2 == 1) {
        //         //     tmpinc = 1;
        //         // } else {
        //         //     tmpinc = -1;
        //         // }
        //         tmpinc = (inc_2 == 1) ? 1 : -1;

		// 		// cout << "###" << inc_2 << " " << prec_ext << " " << tmpinc * prec_ext << endl;
		// 		if(tmpinc * prec_ext >= current_err) {
		// 			local_extreme_point[local_extreme_count].x = prec_x;
		// 			local_extreme_point[local_extreme_count].y = prec_ext;
		// 			// if(inc_2 == 1) {
        //             //     local_extreme_point[local_extreme_count].locmm = 1;
        //             // } else {
        //             //     local_extreme_point[local_extreme_count].locmm = -1;
        //             // }
        //             local_extreme_point[local_extreme_count++].locmm = (inc_2 == 1) ? 1 : -1;

		// 			// cout << extreme_point[extreme_count].x << " " << extreme_point[extreme_count].y << " " << extreme_point[extreme_count].locmm << endl;
		// 			// local_extreme_count++;
		// 			// cout << "prec_end : " << prec_end << endl;

        //             // cout << "[" << section_ind << "]: " << local_extreme_count << endl;
		// 		}
		// 		// if(!prec_end) {
        //         //     inc_2 = (-1) * inc_2;
        //         // }
        //         inc_2 *= prec_end ? 1 : -1;
		// 	}
		// 	inc_2 = 0;
		// } else {
			// inc_1 = inc_2;
			// scan_y1 = scan_y2;
			// scan_y2 = chebeval(deg, coeff, scan_2 / chebeval_k) - function_value(scan_2);
			// if(scan_y1 < scan_y2) {
            //     inc_2 = 1;
            // } else if(scan_y1 > scan_y2) {
            //     inc_2 = -1;
            // } else {
            //     inc_2 = 0;
            // }
			// // cout << "inc : " << inc_2 << endl;
			// // if((inc_1 == 1 && inc_2 != 1) || (inc_1 == -1 && inc_2 != -1) || inc_1 == 0) {
            // if(inc_1 * inc_2 != 1 || inc_1 == 0) {
			// 	prec_end = false;
			// 	tmp_inc = inc_2;
			// 	prec_x = scan_2;
			// 	while(!prec_end) {
			// 		prec_sc = (prec_x - scan_1) / 2;
			// 		prec_end = true;
			// 		if(inc_1 != 0) {
			// 			// cout << "extreme detected" << endl;
			// 			for(long k = 0; k < 3; k++) {
			// 				detail[k] = scan_1 + (k - 1) * prec_sc;
			// 			}
			// 		} else {
			// 			// cout << "boundary start" << endl;
			// 			for(long k = 0; k < 3; k++) {
			// 				detail[k] = scan_1 + prec_sc * k;
			// 			}
			// 		}
			// 		prec_iter = 0;
			// 		while(prec_iter < params.binary_prec) {
			// 			prec_ext = chebeval(deg, coeff, detail[0] / chebeval_k) - function_value(detail[0]);
			// 			prec_ind = 0;
			// 			for(long k = 1; k < 3; k++) {
			// 				tmp = chebeval(deg, coeff, detail[k] / chebeval_k) - function_value(detail[k]);
			// 				if((inc_1 == 1 && prec_ext < tmp) || (inc_1 == -1 && prec_ext > tmp)) {
			// 					prec_ext = tmp;
			// 					prec_ind = k;
			// 				} else if(inc_1 == 0) {
			// 					if((inc_2 == 1 && prec_ext > tmp) || (inc_2 == -1 && prec_ext < tmp)) {
			// 						prec_ext = tmp;
			// 						prec_ind = k;
			// 					}
			// 				}
			// 			}

			// 			if(inc_1 == 0 && prec_ind != 0) {
            //                 prec_end = false;
            //             }
			// 			prec_x = detail[prec_ind];
			// 			prec_sc = prec_sc / 2;
			// 			// cout << "prec_ind : " << prec_ind << " prec_x : " << prec_x << endl;
			// 			for(long k = 0; k < 3; k++) {
			// 				// if(inc_1 == 0) {
            //                 //     detail[k] = max(prec_x - prec_sc, scan_1) + prec_sc * k;
            //                 // } else {
            //                 //     detail[k] = prec_x - prec_sc + prec_sc * k;
            //                 // }
            //                 detail[k] = (inc_1 == 0) ? max(prec_x - prec_sc, scan_1) + prec_sc * k : prec_x + (k - 1) * prec_sc;
			// 			}
			// 			prec_iter++;
			// 		}
			// 		// if(inc_2 == 1) {
            //         //     tmpinc = -1;
            //         // } else {
            //         //     tmpinc = 1;
            //         // }
            //         tmpinc = (inc_2 == 1) ? -1 : 1;

			// 		// cout << "###" << inc_2 << " " << prec_ext << " " << tmpinc * prec_ext << endl;
			// 		if(tmpinc * prec_ext >= current_err) {
			// 			local_extreme_point[local_extreme_count].x = prec_x;
			// 			local_extreme_point[local_extreme_count].y = prec_ext;
			// 			// if(inc_2 == 1) {
            //             //     local_extreme_point[local_extreme_count].locmm = -1;
            //             // } else {
            //             //     local_extreme_point[local_extreme_count].locmm = +1;
            //             // }
            //             local_extreme_point[local_extreme_count++].locmm = (inc_2 == 1) ? -1 : 1;

			// 			// cout << extreme_point[local_extreme_count].x << " " << extreme_point[local_extreme_count].y << " " << extreme_point[local_extreme_count].locmm << endl;
			// 			// local_extreme_count++;
			// 			// cout << "prec_end : " << prec_end << endl;

            //             // cout << "[" << section_ind << "]: " << local_extreme_count << endl;
			// 		}
			// 		// if(!prec_end) { 
			// 		// 	inc_2 = (-1) * inc_2; 
			// 		// 	// cout << "inc : " << inc_2 << endl;
			// 		// }
            //         inc_2 *= prec_end ? 1 : -1;
			// 	}
			// 	inc_2 = tmp_inc;
			// }
		// }
	}
}

void RemezGeneral::getextreme() {
    Point **local_extreme_point_array = new Point * [section_num];
    for(long i = 0; i < section_num; i++) {
        local_extreme_point_array[i] = new Point[4 * deg];
    }
    long local_extreme_count_array[section_num];

    vector<thread> vec_thr;

    for (long i = 0; i < section_num; i++) {
        vec_thr.emplace_back(&RemezGeneral::getextreme_local, this, local_extreme_point_array[i],
        ref(local_extreme_count_array[i]), i);
        // getextreme_local(local_extreme_point_array[i], local_extreme_count_array[i], i - boundary_K + 1);
    }

    for (auto &t : vec_thr) {
        t.join();
    }

    // cout << "getextreme_local has ended" << endl;

    extreme_count = 0;
    for (long i = 0; i < section_num; i++) {
        for (int j = 0; j < local_extreme_count_array[i]; j++) {
            extreme_point[extreme_count].x = local_extreme_point_array[i][j].x;
            extreme_point[extreme_count].y = local_extreme_point_array[i][j].y;
            extreme_point[extreme_count].locmm = local_extreme_point_array[i][j].locmm;

            extreme_count++;
        }
    }

    // cout << "ok" << endl;
    max_err = 0;
    for (long i = 0; i < extreme_count; i++) {
        // cout << extreme_point[i].x << " " << extreme_point[i].y << " " << extreme_point[i].locmm << endl;
        // if (max_err < abs(extreme_point[i].y)) {
        //     max_err = abs(extreme_point[i].y);
        // }
        max_err = (max_err < abs(extreme_point[i].y)) ? abs(extreme_point[i].y) : max_err;
    }
    sort(extreme_point, extreme_point + extreme_count, xcompare);
    cout << "**" << extreme_count << endl;

    for (long i = 0; i < section_num; i++) {
        delete[] local_extreme_point_array[i];
    }

    delete[] local_extreme_point_array;
    // delete[] local_extreme_count_array;
}

void RemezGeneral::choosemaxs() {
    RR::SetPrecision(params.RR_prec);
    Point *extract = new Point[extreme_count];
    long count = 0, ind = 0;
    max_err = RR(0);
    min_err = RR(1000);

    long *temparray = new long[extreme_count];
    long maxtemp, tempcount = 0;
    RR maxtempvalue;

    while (ind < extreme_count) {
        if (tempcount == 0) {
            temparray[tempcount] = ind;
            tempcount++;
            ind++;
        } else {
            if (extreme_point[ind].locmm * extreme_point[ind - 1].locmm == 1) {
                temparray[tempcount] = ind;
                tempcount++;
                ind++;
            } else {
                maxtempvalue = RR(0);
                for (long i = 0; i < tempcount; i++) {
                    if (maxtempvalue < abs(extreme_point[temparray[i]].y)) {
                        maxtempvalue = abs(extreme_point[temparray[i]].y);
                        maxtemp = temparray[i];
                    }
                }
                extract[count].x = extreme_point[maxtemp].x;
                extract[count].y = extreme_point[maxtemp].y;
                extract[count].locmm = extreme_point[maxtemp].locmm;
                count++;
                tempcount = 0;
            }
        }
    }
    
    maxtempvalue = RR(0);
    for (long i = 0; i < tempcount; i++) {
        if (maxtempvalue < abs(extreme_point[temparray[i]].y)) {
            maxtempvalue = abs(extreme_point[temparray[i]].y);
            maxtemp = temparray[i];
        }
    }
    extract[count].x = extreme_point[maxtemp].x;
    extract[count].y = extreme_point[maxtemp].y;
    extract[count].locmm = extreme_point[maxtemp].locmm;
    count++;
    tempcount = 0;

    cout << "cnt: " << count << endl;
    // cout << count << " " << deg + 2 << endl;
    RR minsum;
    long minindex;
    while (count > deg + 2) {
        minsum = 100000;
        if (count == deg + 3) {
            if (abs(extract[0].y) > abs(extract[count - 1].y)) {
                count--;
            } else {
                for (long i = 0; i < count - 1; i++) {
                    extract[i].x = extract[i + 1].x;
                    extract[i].y = extract[i + 1].y;
                    extract[i].locmm = extract[i + 1].locmm;
                }
                count--;
            }
        } else if (count == deg + 4) {
            for (long i = 0; i < count; i++) {
                if (minsum > abs(extract[i].y) + abs(extract[(i + 1) % count].y)) {
                    minsum = abs(extract[i].y) + abs(extract[(i + 1) % count].y);
                    minindex = i;
                }
            }
            if (minindex == count - 1) {
                for (long i = 0; i < count - 2; i++) {
                    extract[i].x = extract[i + 1].x;
                    extract[i].y = extract[i + 1].y;
                    extract[i].locmm = extract[i + 1].locmm;
                }
                count -= 2;
            } else {
                for (long i = minindex; i < count - 2; i++) {
                    extract[i].x = extract[i + 2].x;
                    extract[i].y = extract[i + 2].y;
                    extract[i].locmm = extract[i + 2].locmm;
                }
                count -= 2;
            }
        } else {
            for (long i = 0; i < count - 1; i++) {
                if (minsum > abs(extract[i].y) + abs(extract[i + 1].y)) {
                    minsum = abs(extract[i].y) + abs(extract[i + 1].y);
                    minindex = i;
                }
            }
            if (minindex == 0) {
                for (long i = 0; i < count - 1; i++) {
                    extract[i].x = extract[i + 1].x;
                    extract[i].y = extract[i + 1].y;
                    extract[i].locmm = extract[i + 1].locmm;
                }
                count--;
            } else if (minindex == count - 2) {
                count--;
            } else {
                for (long i = minindex; i < count - 2; i++) {
                    extract[i].x = extract[i + 2].x;
                    extract[i].y = extract[i + 2].y;
                    extract[i].locmm = extract[i + 2].locmm;
                }
                count -= 2;
            }
        }
    }
    
    for (long i = 0; i < (deg + 2); i++) {
        cout << i << " ";
        sample_point[i].x = extract[i].x;
        sample_point[i].y = function_value(sample_point[i].x);
        // // cout << extract[i].y << endl;
        // if (max_err < abs(extract[i].y)) {
        //     //			cout << maxerr << endl;
        //     max_err = abs(extract[i].y);
        // }
        // if (min_err > abs(extract[i].y)) {
        //     min_err = abs(extract[i].y);
        // }
        max_err = (max_err < abs(extract[i].y)) ? abs(extract[i].y) : max_err;
        min_err = (min_err > abs(extract[i].y)) ? min_err : abs(extract[i].y);
    }

    delete[] extract;
    delete[] temparray;
}

void RemezGeneral::generate_optimal_poly(Polynomial &poly) {

    // better_initialize();
    initialize();

    getcoeffwitherr();
    getextreme();
    // choosemaxs();

    // long it = 0;
    // while ((max_err - min_err) / min_err > approx) {
    //     getcoeffwitherr();
    //     getextreme();
    //     choosemaxs();
    //     it++;
    //     // for(long i = 0; i < deg + 2; i++) {
    //     //         cout << sample_point[i].x << endl;
    //     // }
    //     cout << it << "th end" << endl;
    //     cout << max_err << " " << min_err << endl;
    // }

    // // for(long i = 0; i < extreme_count; i++) {
    // //         cout << extreme_point[i].x << " " << extreme_point[i].y << " " << extreme_point[i].locmm << endl;
    // // }
    // // cout << endl;
    
    // poly.set_polynomial(deg, coeff, "cheb");
    // //	showgraph(out, coeff, deg, K, sc);
}

void RemezGeneral::showcoeff() {
    for (int i = 0; i < deg + 1; i++) {
        cout << i << " : " << coeff[i] << endl;
    }
}
