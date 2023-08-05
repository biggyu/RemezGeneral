#include"Polynomial.h"

Polynomial::Polynomial() {}

Polynomial::Polynomial(long _deg) {
	deg = _deg;
	coeff = new RR[deg + 1];
	chebcoeff = new RR[deg + 1];
	for(int i = 0; i < deg + 1; i++) {
		coeff[i] = RR(0);
		chebcoeff[i] = RR(0);
	}
}

Polynomial::Polynomial(long _deg, RR* _coeff, string tag) {
	deg = _deg;
	coeff = new RR[deg + 1];
	chebcoeff = new RR[deg + 1];
	if (tag == "power") {
		for(int i = 0; i < deg + 1; i++) {
			coeff[i] = _coeff[i];
		}
		power_to_cheb();
	}

	else if (tag == "cheb") {
		for(int i = 0; i < deg + 1; i++) {
			chebcoeff[i] = _coeff[i];
		}
		cheb_to_power();
	}
}

Polynomial::~Polynomial() {
	if(coeff) delete[] coeff;
	if(chebcoeff) delete[] chebcoeff;
	if(poly_heap) {
		for (int i = 0; i < heaplen; i++) {
			if(poly_heap[i]) delete poly_heap[i];
		}
	}
}

void Polynomial::set_polynomial(long _deg, RR* _coeff, string tag, double chebval_k) {
	deg = _deg;
	coeff = new RR[deg + 1];
	chebcoeff = new RR[deg + 1];
	if (tag == "power") {
		for(int i = 0; i < deg + 1; i++) {
			coeff[i] = _coeff[i];
		}
		power_to_cheb();
	}

	else if (tag == "cheb") {
		for(int i = 0; i < deg + 1; i++) {
			chebcoeff[i] = _coeff[i];
		}
		cheb_to_power(chebval_k);
	}
}

void Polynomial::set_polynomial(long _deg, RR* _coeff, string tag) {
	deg = _deg;
	coeff = new RR[deg + 1];
	chebcoeff = new RR[deg + 1];
	if (tag == "power") {
		for(int i = 0; i < deg + 1; i++) {
			coeff[i] = _coeff[i];
		}
		power_to_cheb();
	}

	else if (tag == "cheb") {
		for(int i = 0; i < deg + 1; i++) {
			chebcoeff[i] = _coeff[i];
		}
		cheb_to_power();
	}
}

void Polynomial::set_zero_polynomial(long _deg) {
	deg = _deg;
	coeff = new RR[deg + 1]{};
	chebcoeff = new RR[deg + 1]{};
}

// void Polynomial::showcoeff() {
// 	for(int i = 0; i < deg + 1; i++) {
// 		cout << "term " << i << " : " << coeff[i] << endl;
// 	}
// 	cout << endl;
// }
void Polynomial::showcoeff() {
    for (int i = 0; i < deg; i++) {
        // cout << i << " : " << coeff[i] << endl;
		cout << coeff[i] << " x^" << i;
		if(coeff[i + 1] > 0) {
			cout << "+ ";
		}
    }
	// cout << coeff[deg] << endl;
    cout << coeff[deg] << " x^" << deg << " " << endl;
}

void Polynomial::showchebcoeff() {
	for(int i = 0; i < deg + 1; i++) {
		cout << "chebterm " << i << " : " << chebcoeff[i] << endl;
	}
	cout << endl;
}
void Polynomial::copy(Polynomial &poly) {
	deg = poly.deg;
	coeff = new RR[deg + 1];
	chebcoeff = new RR[deg + 1];
	for(int i = 0; i < deg + 1; i++) {
		coeff[i] = poly.coeff[i];
		chebcoeff[i] = poly.chebcoeff[i];
	}
}

void Polynomial::power_to_cheb() {
	Polynomial chebbasis, tmp(deg);
	for(int i = 0; i <= deg; i++) {
		tmp.coeff[i] = coeff[i];
	}

	for(int i = deg; i >= 0; i--) {
		chebyshev(chebbasis, i);
		chebcoeff[i] = tmp.coeff[i] / chebbasis.coeff[i];
		for(int j = 0; j <= i; j++) {
			chebbasis.coeff[j] *= chebcoeff[i];
		}
		tmp.subtinplace(chebbasis);
	}
}

void Polynomial::cheb_to_power() {
	Polynomial chebbasis, tmp(deg);
	for(int i = 0; i <= deg; i++) {
		chebyshev(chebbasis, i);
		for(int j = 0; j <= i; j++) {
			chebbasis.coeff[j] *= chebcoeff[i];
		}
		tmp.addinplace(chebbasis);
	}

	for(int i = 0; i <= deg; i++) {
		coeff[i] = tmp.coeff[i];
	}
}
void Polynomial::cheb_to_power(double chebeval_k) {
	Polynomial chebbasis, tmp(deg);
	for(int i = 0; i <= deg; i++) {
		chebyshev(chebbasis, i);
		for(int j = 0; j <= i; j++) {
			chebbasis.coeff[j] *= chebcoeff[i];
		}
		tmp.addinplace(chebbasis);
	}

	for(int i = 0; i <= deg; i++) {
		coeff[i] = tmp.coeff[i] / pow(chebeval_k, i);
	}
}

RR Polynomial::evaluate(RR value) {
	RR rtn = RR(0), term = RR(1);
	for(int i = 0; i <= deg; i++) {
		rtn += coeff[i] * term;
		term *= value;
	}
	
	return rtn;
}

void Polynomial::constmul(RR constant) {
	for (int i = 0; i <= deg; i++) {
		coeff[i] *= constant;
		chebcoeff[i] *= constant;
	}
}

void Polynomial::mulinplace(Polynomial &poly) {
	Polynomial rtn;
	mul(rtn, (*this), poly);
	copy(rtn);	
}

void Polynomial::addinplace(Polynomial &poly) {
	Polynomial rtn;
	add(rtn, (*this), poly);
	copy(rtn);
}

void Polynomial::subtinplace(Polynomial &poly) {
	Polynomial rtn;
	subt(rtn, (*this), poly);
	copy(rtn);	
}

void Polynomial::change_variable_scale(RR scale) {
	RR div = RR(1);
	for (int i = 0; i <= deg; i++) {
		coeff[i] /= div;
		div *= scale;
	}
	power_to_cheb();
}

void Polynomial::generate_poly_heap_manual(long k, long m) {
	heap_k = k;
	heap_m = m;
	heaplen = (1 << (heap_m + 1)) - 1;
	poly_heap = new Polynomial*[heaplen];
	for(long i = 0; i < heaplen; i++) {
		poly_heap[i] = 0;
	}

	poly_heap[0] = new Polynomial();
	poly_heap[0]->copy((*this));

	long chebdeg = heap_k << heap_m;
	long first, last;
	Polynomial chebgiant;
	for(int i = 0; i < heap_m; i++) {
		chebdeg = chebdeg >> 1;
		first = (1 << i) - 1;
		last = (1 << (i + 1)) - 1;
		chebyshev(chebgiant, chebdeg);
		for(int j = first; j < last; j++) {
			if(poly_heap[j]) {
				if(poly_heap[j]->deg < chebdeg) {
					poly_heap[2 * (j + 1)] = new Polynomial();
					poly_heap[2 * (j + 1)]->copy(*poly_heap[j]);
				}

				else {
					poly_heap[2 * (j + 1) - 1] = new Polynomial();
					poly_heap[2 * (j + 1)] = new Polynomial();
					divide_poly(*poly_heap[2 * (j + 1) - 1], *poly_heap[2 * (j + 1)], *poly_heap[j], chebgiant);
					poly_heap[2 * (j + 1) - 1]->power_to_cheb();
					poly_heap[2 * (j + 1)]->power_to_cheb();
				}
			}
		}
	}
}

void Polynomial::generate_poly_heap_odd() {
	oddbabycount(heap_k, heap_m, deg);
	generate_poly_heap_manual(heap_k, heap_m);
}
void Polynomial::generate_poly_heap() {
	babycount(heap_k, heap_m, deg);
	generate_poly_heap_manual(heap_k, heap_m);
}

void Polynomial::write_heap_to_file(ofstream &out) {
	out << heaplen << endl;

	for (int index = 0; index < heaplen; index++) {
		if(poly_heap[index]) {
			out << index << " " << poly_heap[index]->deg << endl;
			for (int i = 0; i <= poly_heap[index]->deg; i++) {
				out << poly_heap[index]->chebcoeff[i] << endl;
			}
			out << endl;
		}
	}
}

void Polynomial::read_heap_from_file(ifstream &in) {
	long index = 0, in_deg;

	in >> heaplen;

	poly_heap = new Polynomial*[heaplen];

	for(int i = 0; i < heaplen; i++) {
		poly_heap[i] = 0;
	}

	while(index < heaplen - 1) {
		in >> index >> in_deg;
		poly_heap[index] = new Polynomial(in_deg);
		for(int i = 0; i <= in_deg; i++) {
			in >> poly_heap[index]->chebcoeff[i];
		}
		poly_heap[index]->cheb_to_power();
	}

	copy(*poly_heap[0]);
}

// void Polynomial::homomorphic_poly_evaluation(shared_ptr<SEALContext> &context, CKKSEncoder &encoder, Encryptor &encryptor, Evaluator &evaluator, RelinKeys &relin_keys, double scale, Ciphertext &rtn, Ciphertext &cipher) {
// 	Ciphertext** baby;
// 	Ciphertext addtmp1, addtmp2;

// 	Plaintext tmpplain;
// 	parms_id_type last_parms_id;

//     const auto &modulus = iter(context->first_context_data()->parms().coeff_modulus());
// 	if (deg == 1) {
// 		encoder.encode(to_double(coeff[1]), cipher.scale(), tmpplain);
//         if (!tmpplain.is_zero()) {
//             evaluator.mod_switch_to_inplace(tmpplain, cipher.parms_id());
//             evaluator.multiply_plain(cipher, tmpplain, rtn);
//             evaluator.rescale_to_next_inplace(rtn);
//         }
//         else encryptor.encrypt_zero(rtn);

// 		encoder.encode(to_double(coeff[0]), rtn.scale(), tmpplain);
// 		last_parms_id = rtn.parms_id();
// 		evaluator.mod_switch_to_inplace(tmpplain, last_parms_id); 
// 		evaluator.add_plain_inplace(rtn, tmpplain);
// 	}

// 	else if(deg == 2) {
// 		Ciphertext squared_cipher;
// 		evaluator.square(cipher, squared_cipher);
// 		evaluator.relinearize_inplace(squared_cipher, relin_keys);
// 		evaluator.rescale_to_next_inplace(squared_cipher);

// 		encoder.encode(to_double(coeff[2]), squared_cipher.scale(), tmpplain);
//         if (!tmpplain.is_zero()) {
//             evaluator.mod_switch_to_inplace(tmpplain, squared_cipher.parms_id());
//             evaluator.multiply_plain_inplace(squared_cipher, tmpplain);
//             evaluator.rescale_to_next_inplace(squared_cipher);
//         }

// 		encoder.encode(to_double(coeff[1]), cipher.scale(), tmpplain);

//         if (!tmpplain.is_zero()) {
//             evaluator.mod_switch_to_inplace(tmpplain, cipher.parms_id());
//             evaluator.multiply_plain(cipher, tmpplain, rtn);
//             evaluator.rescale_to_next_inplace(rtn);
//         }

// 		encoder.encode(to_double(coeff[0]), rtn.scale(), tmpplain);
// 		last_parms_id = rtn.parms_id();
// 		evaluator.mod_switch_to_inplace(tmpplain, last_parms_id);
// 		evaluator.add_plain_inplace(rtn, tmpplain);
		
// 		make_modulus_equal(context, encoder, evaluator, scale, rtn, squared_cipher, addtmp1, addtmp2);
// 		evaluator.add(addtmp1, addtmp2, rtn);
// 	}

// 	else if(deg == 3) {
// 		Ciphertext squared_cipher, cubic_cipher;
// 		rtn = cipher;

// 		evaluator.square(cipher, squared_cipher);
// 		evaluator.relinearize_inplace(squared_cipher, relin_keys);
// 		evaluator.rescale_to_next_inplace(squared_cipher);

// 		encoder.encode(to_double(coeff[3]), cipher.scale(), tmpplain);
// 		evaluator.mod_switch_to_inplace(tmpplain, cipher.parms_id());
// 		evaluator.multiply_plain(cipher, tmpplain, cubic_cipher);
// 		evaluator.rescale_to_next_inplace(cubic_cipher);

// 		make_modulus_equal(context, encoder, evaluator, scale, cubic_cipher, squared_cipher, addtmp1, addtmp2);
// 		evaluator.multiply(addtmp1, addtmp2, cubic_cipher);
// 		evaluator.relinearize_inplace(cubic_cipher, relin_keys);
// 		evaluator.rescale_to_next_inplace(cubic_cipher);

// 		encoder.encode(to_double(coeff[2]), squared_cipher.scale(), tmpplain);
// 		evaluator.mod_switch_to_inplace(tmpplain, squared_cipher.parms_id());
// 		evaluator.multiply_plain_inplace(squared_cipher, tmpplain);
// 		evaluator.rescale_to_next_inplace(squared_cipher);

// 		encoder.encode(to_double(coeff[1]), cipher.scale(), tmpplain);
// 		evaluator.mod_switch_to_inplace(tmpplain, cipher.parms_id());
// 		evaluator.multiply_plain(cipher, tmpplain, rtn);
// 		evaluator.rescale_to_next_inplace(rtn);

// 		encoder.encode(to_double(coeff[0]), rtn.scale(), tmpplain);
// 		last_parms_id = rtn.parms_id();
// 		evaluator.mod_switch_to_inplace(tmpplain, last_parms_id);
// 		evaluator.add_plain_inplace(rtn, tmpplain);
		
// 		make_modulus_equal(context, encoder, evaluator, scale, rtn, squared_cipher, addtmp1, addtmp2);
// 		evaluator.add(addtmp1, addtmp2, rtn);

// 		make_modulus_equal(context, encoder, evaluator, scale, rtn, cubic_cipher, addtmp1, addtmp2);
// 		evaluator.add(addtmp1, addtmp2, rtn);
// 	}
// 	else {
// 		baby = new Ciphertext*[heap_k];
// 		for (int i = 0; i < heap_k; i++) {
// 			baby[i] = 0;
// 		}
// 		baby[1] = new Ciphertext();
// 		*baby[1] = cipher;

// 		Polynomial tmppoly;

// 		for (int i = 1; i <= floor(log(heap_k - 1) / log(2)); i++) {
// 			baby[(1 << i)] = new Ciphertext();
// 			evaluator.square(*baby[(1 << (i - 1))], *baby[(1 << i)]);
// 			evaluator.relinearize_inplace(*baby[(1 << i)], relin_keys);
// 			evaluator.rescale_to_next_inplace(*baby[(1 << i)]);

// 			evaluator.add_inplace(*baby[(1 << i)], *baby[(1 << i)]);
// 			encoder.encode(-1.0, baby[(1 << i)]->scale(), tmpplain);
// 			evaluator.mod_switch_to_inplace(tmpplain, baby[(1 << i)]->parms_id());
// 			evaluator.add_plain_inplace(*baby[(1 << i)], tmpplain);

// 			chebyshev(tmppoly, (1 << i));
// 		}
// 		long lpow2, res, diff;
// 		Ciphertext tmp;
// 		for (int i = 1; i < heap_k; i++) {
// 			if(!baby[i]) {
// 				lpow2 = (1 << (int)floor(log(i) / log(2)));
// 				res = i - lpow2;
// 				diff = abs(lpow2 - res);

// 				baby[i] = new Ciphertext();
// 				make_modulus_equal(context, encoder, evaluator, scale, *baby[lpow2], *baby[res], addtmp1, addtmp2);

// 				evaluator.multiply(addtmp1, addtmp2, *baby[i]);
// 				evaluator.relinearize_inplace(*baby[i], relin_keys);
// 				evaluator.rescale_to_next_inplace(*baby[i]);

// 				evaluator.add_inplace(*baby[i], *baby[i]);
				
// 				make_modulus_equal(context, encoder, evaluator, scale, *baby[diff], *baby[i], addtmp1, addtmp2);
// 				evaluator.sub(addtmp2, addtmp1, *baby[i]);

// 				chebyshev(tmppoly, i);
// 			}
// 		}

// 		Ciphertext** giant = new Ciphertext*[heap_m];
// 		giant[0] = new Ciphertext();
// 		lpow2 = (1 << ((int)ceil(log(heap_k) / log(2)) - 1));
// 		res = heap_k - lpow2;
// 		diff = abs(lpow2 - res);

// 		if (res == 0) {
// 			*giant[0] = *baby[lpow2];
// 		}

// 		else if (diff == 0) {
// 			evaluator.square(*baby[lpow2], *giant[0]);
// 			evaluator.relinearize_inplace(*giant[0], relin_keys);
// 			evaluator.rescale_to_next_inplace(*giant[0]);

// 			evaluator.add_inplace(*giant[0], *giant[0]);
// 			encoder.encode(-1.0, giant[0]->scale(), tmpplain);
// 			last_parms_id = giant[0]->parms_id();
// 			evaluator.mod_switch_to_inplace(tmpplain, last_parms_id);
// 			evaluator.add_plain_inplace(*giant[0], tmpplain);
// 		}
		
// 		else {
// 			make_modulus_equal(context, encoder, evaluator, scale, *baby[lpow2], *baby[res], addtmp1, addtmp2);
// 			evaluator.multiply(addtmp1, addtmp2, *giant[0]);
// 			evaluator.relinearize_inplace(*giant[0], relin_keys);
// 			evaluator.rescale_to_next_inplace(*giant[0]);
// 			evaluator.add_inplace(*giant[0], *giant[0]);
// 			make_modulus_equal(context, encoder, evaluator, scale, *baby[diff], *giant[0], addtmp1, addtmp2);
// 			evaluator.sub(addtmp2, addtmp1, *giant[0]);
// 		}

// 		for (int i = 1; i < heap_m; i++) {
// 			giant[i] = new Ciphertext();
// 			evaluator.square(*giant[i - 1], *giant[i]);
// 			evaluator.relinearize_inplace(*giant[i], relin_keys);
// 			evaluator.rescale_to_next_inplace(*giant[i]);
// 			evaluator.add_inplace(*giant[i], *giant[i]);
			
// 			encoder.encode(-1.0, giant[i]->scale(), tmpplain);
// 			last_parms_id = giant[i]->parms_id();
// 			evaluator.mod_switch_to_inplace(tmpplain, last_parms_id);
// 			evaluator.add_plain_inplace(*giant[i], tmpplain);
// 		}

//         size_t* mod_index_heap = new size_t[(1 << (heap_m + 1)) - 1]{};

//         size_t initial_mod_index = context->get_context_data(cipher.parms_id())->chain_index();

// 		long heapfirst = (1 << heap_m) - 1;
// 		long heaplast = (1 << (heap_m + 1)) - 1;

//         for (int i = heapfirst; i < heaplast; i++) {
//             if (poly_heap[i]) {
//                 if (i % 2 == 1) mod_index_heap[i] = max((size_t)ceil(log2(poly_heap[i]->deg)) + 1, (size_t)ceil(log2(poly_heap[i + 1]->deg + 1))); 
//                 else mod_index_heap[i] = (size_t)ceil(log2(poly_heap[i]->deg)) + 1; 
//             }
//         }

// 		long depth = heap_m;

// 		while (depth != 1) {
// 			depth--;
// 			heapfirst = (1 << depth) - 1;
// 			heaplast = (1 << (depth + 1)) - 1;
// 			for (int i = heapfirst; i < heaplast; i++) {		
//                 if (poly_heap[i]) {
//                     if (i % 2 == 1) {
//                         if (poly_heap[2 * i + 1]) mod_index_heap[i] = max(mod_index_heap[2 * i + 1] + 1, mod_index_heap[2 * i + 2]);
//                         else mod_index_heap[i] = mod_index_heap[2 * i + 2];
//                         mod_index_heap[i] = max(mod_index_heap[i], (size_t)ceil(log2(poly_heap[i + 1]->deg + 1)));
//                     }

//                     else {
//                         mod_index_heap[i] = max(mod_index_heap[2 * i + 1] + 1, mod_index_heap[2 * i + 2]);
//                     }
//                 }
// 			}
// 		}

//         double* res_err_heap = new double[(1 << (heap_m + 1)) - 1]{};
        
// 		depth = heap_m + 1;
// 		long gindex = 0;

// 		while (depth != 1) {
// 			depth--;
// 			heapfirst = (1 << depth) - 1;
// 			heaplast = (1 << (depth + 1)) - 1;
// 			for (int i = heapfirst; i < heaplast; i++) {		
//                 if (poly_heap[i]) {
//                     mod_index_heap[i] = initial_mod_index - mod_index_heap[i];
//                     if (i % 2 == 1) res_err_heap[i] = giant[gindex]->scale() / (double)modulus[mod_index_heap[i]].value(); 
//                     else res_err_heap[i] = 1.0;
//                     // cout << "res_err " << i << " : " << res_err_heap[i] << endl;
//                 }
// 			}
//             gindex++;
// 		}

// 		Ciphertext** cipherheap = new Ciphertext*[(1 << (heap_m + 1)) - 1];
// 		for (int i = 0; i < (1 << (heap_m + 1)) - 1; i++) {
// 			cipherheap[i] = 0;
// 		}

// 		heapfirst = (1 << heap_m) - 1;
// 		heaplast = (1 << (heap_m + 1)) - 1;
// 		for (int i = heapfirst; i < heaplast; i++) {
// 			if(poly_heap[i]) {
// 				cipherheap[i] = new Ciphertext();
				
// 				encryptor.encrypt_zero(*cipherheap[i]);

// 				encoder.encode(to_double(poly_heap[i]->chebcoeff[0]), scale, tmpplain);
// 				cipherheap[i]->scale() = scale;
// 				evaluator.add_plain_inplace(*cipherheap[i], tmpplain);

// 				for (int j = 1; j <= poly_heap[i]->deg; j++) {
// 					encoder.encode(to_double(poly_heap[i]->chebcoeff[j]), baby[j]->scale(), tmpplain);
//                     if (!tmpplain.is_zero()) {
//                         if (j < heap_k) {
//                             last_parms_id = baby[j]->parms_id();
//                             evaluator.mod_switch_to_inplace(tmpplain, last_parms_id);
//                             evaluator.multiply_plain(*baby[j], tmpplain, tmp);
//                         }

//                         else {
//                             encoder.encode(to_double(poly_heap[i]->chebcoeff[j]), giant[0]->scale(), tmpplain);
//                             last_parms_id = giant[0]->parms_id();
//                             evaluator.mod_switch_to_inplace(tmpplain, last_parms_id);
//                             evaluator.multiply_plain(*giant[0], tmpplain, tmp);
//                         }

//                         evaluator.relinearize_inplace(tmp, relin_keys);
//                         evaluator.rescale_to_next_inplace(tmp);
//                         make_modulus_equal(context, encoder, evaluator, scale, *cipherheap[i], tmp, addtmp1, addtmp2);
//                         evaluator.add(addtmp1, addtmp2, *cipherheap[i]);
//                     }
// 				}
// 			}
// 		}
// 		depth = heap_m;
// 		gindex = 0;

// 		while (depth != 0) {
// 			depth--;
// 			heapfirst = (1 << depth) - 1;
// 			heaplast = (1 << (depth + 1)) - 1;
// 			for (int i = heapfirst; i < heaplast; i++) {		
// 				if(poly_heap[i]) {
// 					cipherheap[i] = new Ciphertext();
// 					if (!cipherheap[2 * (i + 1) - 1]) {
// 						*cipherheap[i] = *cipherheap[2 * (i + 1)];
// 					}
					
// 					else {
// 						make_modulus_equal(context, encoder, evaluator, scale, *cipherheap[2 * (i + 1) - 1], *giant[gindex], addtmp1, addtmp2);
// 						evaluator.multiply(addtmp1, addtmp2, *cipherheap[i]);
// 						evaluator.relinearize_inplace(*cipherheap[i], relin_keys);
// 						evaluator.rescale_to_next_inplace(*cipherheap[i]);
// 						make_modulus_equal(context, encoder, evaluator, scale, *cipherheap[i], *cipherheap[2 * (i + 1)], addtmp1, addtmp2);
// 						evaluator.add(addtmp1, addtmp2, *cipherheap[i]);
// 					}

// 				}

// 			}
// 			gindex++;
// 		}

// 		rtn = *cipherheap[0];
// 		for (int i = 0; i < heap_k; i++) {
// 			delete baby[i];
// 		}
		
// 		delete[] baby;

// 		for (int i = 0; i < heap_m; i++) {
// 			delete giant[i];
// 		}

// 		delete[] giant;

// 		for (int i = 0; i < (1 << (heap_m + 1)) - 1; i++) {
// 			delete cipherheap[i];
// 		}

// 		delete[] cipherheap;
// 	}
// }

void mul(Polynomial &rtn, Polynomial &a, Polynomial &b) {
	rtn.set_zero_polynomial(a.deg + b.deg);
	for(long i = 0; i <= rtn.deg; i++) {
		for(long j = 0; j <= i; j++) {
			if(j <= a.deg && i - j <= b.deg) {
				rtn.coeff[i] += a.coeff[j] * b.coeff[i - j];
			}
		}
	}
}

void add(Polynomial &rtn, Polynomial &a, Polynomial &b) {
	if(a.deg >= b.deg) {
		rtn.copy(a);
		for(int i = 0; i <= b.deg; i++) {
			rtn.coeff[i] += b.coeff[i];
		}
	}

	else {
		rtn.copy(b);
		for(int i = 0; i < a.deg; i++) {
			rtn.coeff[i] += a.coeff[i];
		}
	}
}

void subt(Polynomial &rtn, Polynomial &a, Polynomial &b) {
	if(a.deg >= b.deg) {
		rtn.copy(a);
		for(int i = 0; i <= b.deg; i++) {
			rtn.coeff[i] -= b.coeff[i];
		}
	}

	else {
		rtn.copy(b);
		for(int i = 0; i < b.deg; i++) {
			rtn.coeff[i] *= RR(-1);
		}

		for(int i = 0; i < a.deg; i++) {
			rtn.coeff[i] += a.coeff[i];
		}
	}
}


void divide_poly(Polynomial &quotient, Polynomial &remainder, Polynomial &target, Polynomial &divider) {
	if(target.deg < divider.deg) {
		quotient.set_zero_polynomial(0);
		quotient.coeff[0] = 0;

		remainder.copy(target);
	}

	else {
		quotient.set_zero_polynomial(target.deg - divider.deg);
		
		Polynomial tmp(target.deg, target.coeff, "power");
		RR currcoeff;
		long currdeg = target.deg;
		for(long i = quotient.deg; i >= 0; i--) {
			currcoeff = tmp.coeff[currdeg] / divider.coeff[divider.deg];
			quotient.coeff[i] = currcoeff;
			tmp.coeff[currdeg] = 0;
			for(int j = 0; j < divider.deg; j++) {
				tmp.coeff[j + i] -= divider.coeff[j] * currcoeff;
			}
			currdeg--;
			
		}
		remainder.set_polynomial(divider.deg - 1, tmp.coeff, "power");
	}
}

void chebyshev(Polynomial &rtn, long deg) {
	if(deg == 0) {
		rtn.set_zero_polynomial(0);
		rtn.coeff[0] = RR(1);
	}

	else if(deg == 1) {
		rtn.set_zero_polynomial(1);
		rtn.coeff[1] = RR(1);
	}

	else {
		Polynomial iden2(1);
		iden2.coeff[1] = RR(2);

		Polynomial tmp1(0), tmp2(1), tmp3;

		tmp1.coeff[0] = RR(1);
		tmp2.coeff[1] = RR(1);

		for(int i = 2; i <= deg; i++) {
			mul(tmp3, iden2, tmp2);
			subt(rtn, tmp3, tmp1);
			tmp1.copy(tmp2);
			tmp2.copy(rtn);
		}
	}
}

void second_chebyshev_times_x_for_sine(Polynomial &rtn, long deg) {
	if(deg == 1) {
		rtn.set_zero_polynomial(1);
		rtn.coeff[1] = RR(1);
	}

	else if(deg == 3) {
		rtn.set_zero_polynomial(3);
		rtn.coeff[1] = RR(3);
		rtn.coeff[3] = RR(-4);
	}

	else if(deg % 2 == 1 && deg > 0){
		Polynomial mul_fac(2);
		mul_fac.coeff[0] = RR(2);
		mul_fac.coeff[2] = RR(-4);

		Polynomial iden(1);
		iden.coeff[1] = RR(1);

		Polynomial tmp1(0), tmp2(2), tmp3;

		tmp1.coeff[0] = RR(1);
		tmp2.coeff[0] = RR(3);
		tmp2.coeff[2] = RR(-4);

		for(int i = 4; i < deg; i += 2) {
			mul(tmp3, mul_fac, tmp2);
			subt(rtn, tmp3, tmp1);
			tmp1.copy(tmp2);
			tmp2.copy(rtn);
		
		}
		rtn.mulinplace(iden);
	} 
}