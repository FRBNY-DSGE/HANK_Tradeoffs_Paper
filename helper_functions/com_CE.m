
CEs_vec        = CEs(:);
CEs_PROFIT_vec = CEs_PROFIT(:);
CEs_RR_vec     = CEs_RR(:);
CEs_W_vec      = CEs_W(:);
CEs_f_vec      = CEs_f(:);
CEs_GOV_vec    = CEs_GOV(:);
CEs_Q_vec      = CEs_Q(:);

NW_vec = param.q*meshes.a(:)+meshes.b(:);

CE_MEDIAN             = weightedMedian(CEs(:),SS_stats.mu_dist(:)                       /sum(SS_stats.mu_dist(:)));
CE_MEDIAN_Q1          = weightedMedian(CEs(:),SS_stats.mu_dist_Q1(:)                    /sum(SS_stats.mu_dist_Q1(:)));   
CE_MEDIAN_Q2          = weightedMedian(CEs(:),SS_stats.mu_dist_Q2(:)                    /sum(SS_stats.mu_dist_Q2(:)));   
CE_MEDIAN_Q3          = weightedMedian(CEs(:),SS_stats.mu_dist_Q3(:)                    /sum(SS_stats.mu_dist_Q3(:)));   
CE_MEDIAN_Q4          = weightedMedian(CEs(:),SS_stats.mu_dist_Q4(:)                    /sum(SS_stats.mu_dist_Q4(:)));   
CE_MEDIAN_Q5          = weightedMedian(CEs(:),SS_stats.mu_dist_Q5(:)                    /sum(SS_stats.mu_dist_Q5(:)));   
CE_MEDIAN_P90         = weightedMedian(CEs(:),SS_stats.mu_dist_P90(:)                   /sum(SS_stats.mu_dist_P90(:)));   
CE_MEDIAN_P99         = weightedMedian(CEs(:),SS_stats.mu_dist_P99(:)                   /sum(SS_stats.mu_dist_P99(:)));   
CE_MEDIAN_P999        = weightedMedian(CEs(:),SS_stats.mu_dist_P999(:)                  /sum(SS_stats.mu_dist_P999(:)));   
CE_MEDIAN_Ent         = weightedMedian(CEs(:),SS_stats.mu_dist_Ent(:)                   /sum(SS_stats.mu_dist_Ent(:)));   
CE_MEDIAN_Emp         = weightedMedian(CEs(:),SS_stats.mu_dist_Emp(:)                   /sum(SS_stats.mu_dist_Emp(:)));   
CE_MEDIAN_Unemp       = weightedMedian(CEs(:),SS_stats.mu_dist_Unemp(:)                 /sum(SS_stats.mu_dist_Unemp(:)));   
CE_MEDIAN_Ent         = weightedMedian(CEs(:),SS_stats.mu_dist_Ent(:).*SS_stats.mu_dist(:)       /sum(SS_stats.mu_dist_Ent(:).*SS_stats.mu_dist(:)));
CE_MEDIAN_Ent_Q1      = weightedMedian(CEs(:),SS_stats.mu_dist_Ent(:).*SS_stats.mu_dist_Q1(:)    /sum(SS_stats.mu_dist_Ent(:).*SS_stats.mu_dist_Q1(:)));   
CE_MEDIAN_Ent_Q2      = weightedMedian(CEs(:),SS_stats.mu_dist_Ent(:).*SS_stats.mu_dist_Q2(:)    /sum(SS_stats.mu_dist_Ent(:).*SS_stats.mu_dist_Q2(:)));   
CE_MEDIAN_Ent_Q3      = weightedMedian(CEs(:),SS_stats.mu_dist_Ent(:).*SS_stats.mu_dist_Q3(:)    /sum(SS_stats.mu_dist_Ent(:).*SS_stats.mu_dist_Q3(:)));   
CE_MEDIAN_Ent_Q4      = weightedMedian(CEs(:),SS_stats.mu_dist_Ent(:).*SS_stats.mu_dist_Q4(:)    /sum(SS_stats.mu_dist_Ent(:).*SS_stats.mu_dist_Q4(:)));   
CE_MEDIAN_Ent_Q5      = weightedMedian(CEs(:),SS_stats.mu_dist_Ent(:).*SS_stats.mu_dist_Q5(:)    /sum(SS_stats.mu_dist_Ent(:).*SS_stats.mu_dist_Q5(:)));   
CE_MEDIAN_Ent_P90     = weightedMedian(CEs(:),SS_stats.mu_dist_Ent(:).*SS_stats.mu_dist_P90(:)   /sum(SS_stats.mu_dist_Ent(:).*SS_stats.mu_dist_P90(:)));   
CE_MEDIAN_Ent_P99     = weightedMedian(CEs(:),SS_stats.mu_dist_Ent(:).*SS_stats.mu_dist_P99(:)   /sum(SS_stats.mu_dist_Ent(:).*SS_stats.mu_dist_P99(:)));   
CE_MEDIAN_Ent_P999    = weightedMedian(CEs(:),SS_stats.mu_dist_Ent(:).*SS_stats.mu_dist_P999(:)  /sum(SS_stats.mu_dist_Ent(:).*SS_stats.mu_dist_P999(:)));   
CE_MEDIAN_Emp         = weightedMedian(CEs(:),SS_stats.mu_dist_Emp(:).*SS_stats.mu_dist(:)       /sum(SS_stats.mu_dist_Emp(:).*SS_stats.mu_dist(:)));
CE_MEDIAN_Emp_Q1      = weightedMedian(CEs(:),SS_stats.mu_dist_Emp(:).*SS_stats.mu_dist_Q1(:)    /sum(SS_stats.mu_dist_Emp(:).*SS_stats.mu_dist_Q1(:)));   
CE_MEDIAN_Emp_Q2      = weightedMedian(CEs(:),SS_stats.mu_dist_Emp(:).*SS_stats.mu_dist_Q2(:)    /sum(SS_stats.mu_dist_Emp(:).*SS_stats.mu_dist_Q2(:)));   
CE_MEDIAN_Emp_Q3      = weightedMedian(CEs(:),SS_stats.mu_dist_Emp(:).*SS_stats.mu_dist_Q3(:)    /sum(SS_stats.mu_dist_Emp(:).*SS_stats.mu_dist_Q3(:)));   
CE_MEDIAN_Emp_Q4      = weightedMedian(CEs(:),SS_stats.mu_dist_Emp(:).*SS_stats.mu_dist_Q4(:)    /sum(SS_stats.mu_dist_Emp(:).*SS_stats.mu_dist_Q4(:)));   
CE_MEDIAN_Emp_Q5      = weightedMedian(CEs(:),SS_stats.mu_dist_Emp(:).*SS_stats.mu_dist_Q5(:)    /sum(SS_stats.mu_dist_Emp(:).*SS_stats.mu_dist_Q5(:)));   
CE_MEDIAN_Emp_P90     = weightedMedian(CEs(:),SS_stats.mu_dist_Emp(:).*SS_stats.mu_dist_P90(:)   /sum(SS_stats.mu_dist_Emp(:).*SS_stats.mu_dist_P90(:)));   
CE_MEDIAN_Emp_P99     = weightedMedian(CEs(:),SS_stats.mu_dist_Emp(:).*SS_stats.mu_dist_P99(:)   /sum(SS_stats.mu_dist_Emp(:).*SS_stats.mu_dist_P99(:)));   
CE_MEDIAN_Emp_P999    = weightedMedian(CEs(:),SS_stats.mu_dist_Emp(:).*SS_stats.mu_dist_P999(:)  /sum(SS_stats.mu_dist_Emp(:).*SS_stats.mu_dist_P999(:)));   
CE_MEDIAN_Unemp       = weightedMedian(CEs(:),SS_stats.mu_dist_Unemp(:).*SS_stats.mu_dist(:)     /sum(SS_stats.mu_dist_Unemp(:).*SS_stats.mu_dist(:)));
CE_MEDIAN_Unemp_Q1    = weightedMedian(CEs(:),SS_stats.mu_dist_Unemp(:).*SS_stats.mu_dist_Q1(:)  /sum(SS_stats.mu_dist_Unemp(:).*SS_stats.mu_dist_Q1(:)));   
CE_MEDIAN_Unemp_Q2    = weightedMedian(CEs(:),SS_stats.mu_dist_Unemp(:).*SS_stats.mu_dist_Q2(:)  /sum(SS_stats.mu_dist_Unemp(:).*SS_stats.mu_dist_Q2(:)));   
CE_MEDIAN_Unemp_Q3    = weightedMedian(CEs(:),SS_stats.mu_dist_Unemp(:).*SS_stats.mu_dist_Q3(:)  /sum(SS_stats.mu_dist_Unemp(:).*SS_stats.mu_dist_Q3(:)));   
CE_MEDIAN_Unemp_Q4    = weightedMedian(CEs(:),SS_stats.mu_dist_Unemp(:).*SS_stats.mu_dist_Q4(:)  /sum(SS_stats.mu_dist_Unemp(:).*SS_stats.mu_dist_Q4(:)));   
CE_MEDIAN_Unemp_Q5    = weightedMedian(CEs(:),SS_stats.mu_dist_Unemp(:).*SS_stats.mu_dist_Q5(:)  /sum(SS_stats.mu_dist_Unemp(:).*SS_stats.mu_dist_Q5(:)));   
CE_MEDIAN_Unemp_P90   = weightedMedian(CEs(:),SS_stats.mu_dist_Unemp(:).*SS_stats.mu_dist_P90(:) /sum(SS_stats.mu_dist_Unemp(:).*SS_stats.mu_dist_P90(:)));   
CE_MEDIAN_Unemp_P99   = weightedMedian(CEs(:),SS_stats.mu_dist_Unemp(:).*SS_stats.mu_dist_P99(:) /sum(SS_stats.mu_dist_Unemp(:).*SS_stats.mu_dist_P99(:)));   
CE_MEDIAN_Unemp_P999  = weightedMedian(CEs(:),SS_stats.mu_dist_Unemp(:).*SS_stats.mu_dist_P999(:)/sum(SS_stats.mu_dist_Unemp(:).*SS_stats.mu_dist_P999(:)));   

CE_P10          = sum(CEs_vec.*(NW_vec==w10)      .*mu_dist_Q1(:))       /sum((NW_vec==w10)      .*mu_dist_Q1(:));
CE_P90          = sum(CEs_vec.*(NW_vec==w90)      .*mu_dist_Q5(:))       /sum((NW_vec==w90)      .*mu_dist_Q5(:));
CE_P90toP10_1   = CE_P90/CE_P10;
CE_P90toP10_2   = CE_MEDIAN_Q5/CE_MEDIAN_Q1;

CE_MEDIAN_PROFIT             = sum(CEs_PROFIT_vec.*(CEs_vec==CE_MEDIAN))      /sum((CEs_vec==CE_MEDIAN));
CE_MEDIAN_PROFIT_Q1          = sum(CEs_PROFIT_vec.*(CEs_vec==CE_MEDIAN_Q1))   /sum((CEs_vec==CE_MEDIAN_Q1));
CE_MEDIAN_PROFIT_Q2          = sum(CEs_PROFIT_vec.*(CEs_vec==CE_MEDIAN_Q2))   /sum((CEs_vec==CE_MEDIAN_Q2));
CE_MEDIAN_PROFIT_Q3          = sum(CEs_PROFIT_vec.*(CEs_vec==CE_MEDIAN_Q3))   /sum((CEs_vec==CE_MEDIAN_Q3));
CE_MEDIAN_PROFIT_Q4          = sum(CEs_PROFIT_vec.*(CEs_vec==CE_MEDIAN_Q4))   /sum((CEs_vec==CE_MEDIAN_Q4));
CE_MEDIAN_PROFIT_Q5          = sum(CEs_PROFIT_vec.*(CEs_vec==CE_MEDIAN_Q5))   /sum((CEs_vec==CE_MEDIAN_Q5));
CE_MEDIAN_PROFIT_P90         = sum(CEs_PROFIT_vec.*(CEs_vec==CE_MEDIAN_P90))  /sum((CEs_vec==CE_MEDIAN_P90));
CE_MEDIAN_PROFIT_P99         = sum(CEs_PROFIT_vec.*(CEs_vec==CE_MEDIAN_P99))  /sum((CEs_vec==CE_MEDIAN_P99));
CE_MEDIAN_PROFIT_P999        = sum(CEs_PROFIT_vec.*(CEs_vec==CE_MEDIAN_P999)) /sum((CEs_vec==CE_MEDIAN_P999));
CE_MEDIAN_PROFIT_Ent         = sum(CEs_PROFIT_vec.*(CEs_vec==CE_MEDIAN_Ent))  /sum((CEs_vec==CE_MEDIAN_Ent));
CE_MEDIAN_PROFIT_Emp         = sum(CEs_PROFIT_vec.*(CEs_vec==CE_MEDIAN_Emp))  /sum((CEs_vec==CE_MEDIAN_Emp));
CE_MEDIAN_PROFIT_Unemp       = sum(CEs_PROFIT_vec.*(CEs_vec==CE_MEDIAN_Unemp))/sum((CEs_vec==CE_MEDIAN_Unemp));
CE_MEDIAN_RR                 = sum(CEs_RR_vec    .*(CEs_vec==CE_MEDIAN))      /sum((CEs_vec==CE_MEDIAN));
CE_MEDIAN_RR_Q1              = sum(CEs_RR_vec    .*(CEs_vec==CE_MEDIAN_Q1))   /sum((CEs_vec==CE_MEDIAN_Q1));
CE_MEDIAN_RR_Q2              = sum(CEs_RR_vec    .*(CEs_vec==CE_MEDIAN_Q2))   /sum((CEs_vec==CE_MEDIAN_Q2));
CE_MEDIAN_RR_Q3              = sum(CEs_RR_vec    .*(CEs_vec==CE_MEDIAN_Q3))   /sum((CEs_vec==CE_MEDIAN_Q3));
CE_MEDIAN_RR_Q4              = sum(CEs_RR_vec    .*(CEs_vec==CE_MEDIAN_Q4))   /sum((CEs_vec==CE_MEDIAN_Q4));
CE_MEDIAN_RR_Q5              = sum(CEs_RR_vec    .*(CEs_vec==CE_MEDIAN_Q5))   /sum((CEs_vec==CE_MEDIAN_Q5));
CE_MEDIAN_RR_P90             = sum(CEs_RR_vec    .*(CEs_vec==CE_MEDIAN_P90))  /sum((CEs_vec==CE_MEDIAN_P90));
CE_MEDIAN_RR_P99             = sum(CEs_RR_vec    .*(CEs_vec==CE_MEDIAN_P99))  /sum((CEs_vec==CE_MEDIAN_P99));
CE_MEDIAN_RR_P999            = sum(CEs_RR_vec    .*(CEs_vec==CE_MEDIAN_P999)) /sum((CEs_vec==CE_MEDIAN_P999));
CE_MEDIAN_RR_Ent             = sum(CEs_RR_vec    .*(CEs_vec==CE_MEDIAN_Ent))  /sum((CEs_vec==CE_MEDIAN_Ent));
CE_MEDIAN_RR_Emp             = sum(CEs_RR_vec    .*(CEs_vec==CE_MEDIAN_Emp))  /sum((CEs_vec==CE_MEDIAN_Emp));
CE_MEDIAN_RR_Unemp           = sum(CEs_RR_vec    .*(CEs_vec==CE_MEDIAN_Unemp))/sum((CEs_vec==CE_MEDIAN_Unemp));
CE_MEDIAN_W                  = sum(CEs_W_vec     .*(CEs_vec==CE_MEDIAN))      /sum((CEs_vec==CE_MEDIAN));
CE_MEDIAN_W_Q1               = sum(CEs_W_vec     .*(CEs_vec==CE_MEDIAN_Q1))   /sum((CEs_vec==CE_MEDIAN_Q1));
CE_MEDIAN_W_Q2               = sum(CEs_W_vec     .*(CEs_vec==CE_MEDIAN_Q2))   /sum((CEs_vec==CE_MEDIAN_Q2));
CE_MEDIAN_W_Q3               = sum(CEs_W_vec     .*(CEs_vec==CE_MEDIAN_Q3))   /sum((CEs_vec==CE_MEDIAN_Q3));
CE_MEDIAN_W_Q4               = sum(CEs_W_vec     .*(CEs_vec==CE_MEDIAN_Q4))   /sum((CEs_vec==CE_MEDIAN_Q4));
CE_MEDIAN_W_Q5               = sum(CEs_W_vec     .*(CEs_vec==CE_MEDIAN_Q5))   /sum((CEs_vec==CE_MEDIAN_Q5));
CE_MEDIAN_W_P90              = sum(CEs_W_vec     .*(CEs_vec==CE_MEDIAN_P90))  /sum((CEs_vec==CE_MEDIAN_P90));
CE_MEDIAN_W_P99              = sum(CEs_W_vec     .*(CEs_vec==CE_MEDIAN_P99))  /sum((CEs_vec==CE_MEDIAN_P99));
CE_MEDIAN_W_P999             = sum(CEs_W_vec     .*(CEs_vec==CE_MEDIAN_P999)) /sum((CEs_vec==CE_MEDIAN_P999));
CE_MEDIAN_W_Ent              = sum(CEs_W_vec     .*(CEs_vec==CE_MEDIAN_Ent))  /sum((CEs_vec==CE_MEDIAN_Ent));
CE_MEDIAN_W_Emp              = sum(CEs_W_vec     .*(CEs_vec==CE_MEDIAN_Emp))  /sum((CEs_vec==CE_MEDIAN_Emp));
CE_MEDIAN_W_Unemp            = sum(CEs_W_vec     .*(CEs_vec==CE_MEDIAN_Unemp))/sum((CEs_vec==CE_MEDIAN_Unemp));
CE_MEDIAN_f                  = sum(CEs_f_vec     .*(CEs_vec==CE_MEDIAN))      /sum((CEs_vec==CE_MEDIAN));
CE_MEDIAN_f_Q1               = sum(CEs_f_vec     .*(CEs_vec==CE_MEDIAN_Q1))   /sum((CEs_vec==CE_MEDIAN_Q1));
CE_MEDIAN_f_Q2               = sum(CEs_f_vec     .*(CEs_vec==CE_MEDIAN_Q2))   /sum((CEs_vec==CE_MEDIAN_Q2));
CE_MEDIAN_f_Q3               = sum(CEs_f_vec     .*(CEs_vec==CE_MEDIAN_Q3))   /sum((CEs_vec==CE_MEDIAN_Q3));
CE_MEDIAN_f_Q4               = sum(CEs_f_vec     .*(CEs_vec==CE_MEDIAN_Q4))   /sum((CEs_vec==CE_MEDIAN_Q4));
CE_MEDIAN_f_Q5               = sum(CEs_f_vec     .*(CEs_vec==CE_MEDIAN_Q5))   /sum((CEs_vec==CE_MEDIAN_Q5));
CE_MEDIAN_f_P90              = sum(CEs_f_vec     .*(CEs_vec==CE_MEDIAN_P90))  /sum((CEs_vec==CE_MEDIAN_P90));
CE_MEDIAN_f_P99              = sum(CEs_f_vec     .*(CEs_vec==CE_MEDIAN_P99))  /sum((CEs_vec==CE_MEDIAN_P99));
CE_MEDIAN_f_P999             = sum(CEs_f_vec     .*(CEs_vec==CE_MEDIAN_P999)) /sum((CEs_vec==CE_MEDIAN_P999));
CE_MEDIAN_f_Ent              = sum(CEs_f_vec     .*(CEs_vec==CE_MEDIAN_Ent))  /sum((CEs_vec==CE_MEDIAN_Ent));
CE_MEDIAN_f_Emp              = sum(CEs_f_vec     .*(CEs_vec==CE_MEDIAN_Emp))  /sum((CEs_vec==CE_MEDIAN_Emp));
CE_MEDIAN_f_Unemp            = sum(CEs_f_vec     .*(CEs_vec==CE_MEDIAN_Unemp))/sum((CEs_vec==CE_MEDIAN_Unemp));
CE_MEDIAN_GOV                = sum(CEs_GOV_vec   .*(CEs_vec==CE_MEDIAN))      /sum((CEs_vec==CE_MEDIAN));
CE_MEDIAN_GOV_Q1             = sum(CEs_GOV_vec   .*(CEs_vec==CE_MEDIAN_Q1))   /sum((CEs_vec==CE_MEDIAN_Q1));
CE_MEDIAN_GOV_Q2             = sum(CEs_GOV_vec   .*(CEs_vec==CE_MEDIAN_Q2))   /sum((CEs_vec==CE_MEDIAN_Q2));
CE_MEDIAN_GOV_Q3             = sum(CEs_GOV_vec   .*(CEs_vec==CE_MEDIAN_Q3))   /sum((CEs_vec==CE_MEDIAN_Q3));
CE_MEDIAN_GOV_Q4             = sum(CEs_GOV_vec   .*(CEs_vec==CE_MEDIAN_Q4))   /sum((CEs_vec==CE_MEDIAN_Q4));
CE_MEDIAN_GOV_Q5             = sum(CEs_GOV_vec   .*(CEs_vec==CE_MEDIAN_Q5))   /sum((CEs_vec==CE_MEDIAN_Q5));
CE_MEDIAN_GOV_P90            = sum(CEs_GOV_vec   .*(CEs_vec==CE_MEDIAN_P90))  /sum((CEs_vec==CE_MEDIAN_P90));
CE_MEDIAN_GOV_P99            = sum(CEs_GOV_vec   .*(CEs_vec==CE_MEDIAN_P99))  /sum((CEs_vec==CE_MEDIAN_P99));
CE_MEDIAN_GOV_P999           = sum(CEs_GOV_vec   .*(CEs_vec==CE_MEDIAN_P999)) /sum((CEs_vec==CE_MEDIAN_P999));
CE_MEDIAN_GOV_Ent            = sum(CEs_GOV_vec   .*(CEs_vec==CE_MEDIAN_Ent))  /sum((CEs_vec==CE_MEDIAN_Ent));
CE_MEDIAN_GOV_Emp            = sum(CEs_GOV_vec   .*(CEs_vec==CE_MEDIAN_Emp))  /sum((CEs_vec==CE_MEDIAN_Emp));
CE_MEDIAN_GOV_Unemp          = sum(CEs_GOV_vec   .*(CEs_vec==CE_MEDIAN_Unemp))/sum((CEs_vec==CE_MEDIAN_Unemp));
CE_MEDIAN_Q                  = sum(CEs_Q_vec     .*(CEs_vec==CE_MEDIAN))      /sum((CEs_vec==CE_MEDIAN));
CE_MEDIAN_Q_Q1               = sum(CEs_Q_vec     .*(CEs_vec==CE_MEDIAN_Q1))   /sum((CEs_vec==CE_MEDIAN_Q1));
CE_MEDIAN_Q_Q2               = sum(CEs_Q_vec     .*(CEs_vec==CE_MEDIAN_Q2))   /sum((CEs_vec==CE_MEDIAN_Q2));
CE_MEDIAN_Q_Q3               = sum(CEs_Q_vec     .*(CEs_vec==CE_MEDIAN_Q3))   /sum((CEs_vec==CE_MEDIAN_Q3));
CE_MEDIAN_Q_Q4               = sum(CEs_Q_vec     .*(CEs_vec==CE_MEDIAN_Q4))   /sum((CEs_vec==CE_MEDIAN_Q4));
CE_MEDIAN_Q_Q5               = sum(CEs_Q_vec     .*(CEs_vec==CE_MEDIAN_Q5))   /sum((CEs_vec==CE_MEDIAN_Q5));
CE_MEDIAN_Q_P90              = sum(CEs_Q_vec     .*(CEs_vec==CE_MEDIAN_P90))  /sum((CEs_vec==CE_MEDIAN_P90));
CE_MEDIAN_Q_P99              = sum(CEs_Q_vec     .*(CEs_vec==CE_MEDIAN_P99))  /sum((CEs_vec==CE_MEDIAN_P99));
CE_MEDIAN_Q_P999             = sum(CEs_Q_vec     .*(CEs_vec==CE_MEDIAN_P999)) /sum((CEs_vec==CE_MEDIAN_P999));
CE_MEDIAN_Q_Ent              = sum(CEs_Q_vec     .*(CEs_vec==CE_MEDIAN_Ent))  /sum((CEs_vec==CE_MEDIAN_Ent));
CE_MEDIAN_Q_Emp              = sum(CEs_Q_vec     .*(CEs_vec==CE_MEDIAN_Emp))  /sum((CEs_vec==CE_MEDIAN_Emp));
CE_MEDIAN_Q_Unemp            = sum(CEs_Q_vec     .*(CEs_vec==CE_MEDIAN_Unemp))/sum((CEs_vec==CE_MEDIAN_Unemp));


CE       = sum( SS_stats.mu_dist(:)      .*CEs(:) )/sum( SS_stats.mu_dist(:) );
CE_B01   = sum( SS_stats.mu_dist_B01(:)  .*CEs(:) )/sum( SS_stats.mu_dist_B01(:) );
CE_B1    = sum( SS_stats.mu_dist_B1(:)   .*CEs(:) )/sum( SS_stats.mu_dist_B1(:) );
CE_B10   = sum( SS_stats.mu_dist_B10(:)  .*CEs(:) )/sum( SS_stats.mu_dist_B10(:) );
CE_Q1    = sum( SS_stats.mu_dist_Q1(:)   .*CEs(:) )/sum( SS_stats.mu_dist_Q1(:) );
CE_Q2    = sum( SS_stats.mu_dist_Q2(:)   .*CEs(:) )/sum( SS_stats.mu_dist_Q2(:) );
CE_Q3    = sum( SS_stats.mu_dist_Q3(:)   .*CEs(:) )/sum( SS_stats.mu_dist_Q3(:) );
CE_Q4    = sum( SS_stats.mu_dist_Q4(:)   .*CEs(:) )/sum( SS_stats.mu_dist_Q4(:) );
CE_Q5    = sum( SS_stats.mu_dist_Q5(:)   .*CEs(:) )/sum( SS_stats.mu_dist_Q5(:) );
CE_P90   = sum( SS_stats.mu_dist_P90(:)  .*CEs(:) )/sum( SS_stats.mu_dist_P90(:) );
CE_P99   = sum( SS_stats.mu_dist_P99(:)  .*CEs(:) )/sum( SS_stats.mu_dist_P99(:) );
CE_P999  = sum( SS_stats.mu_dist_P999(:) .*CEs(:) )/sum( SS_stats.mu_dist_P999(:) );
CE_Ent   = sum( SS_stats.mu_dist_Ent(:)  .*CEs(:) )/sum( SS_stats.mu_dist_Ent(:) );
CE_Emp   = sum( SS_stats.mu_dist_Emp(:)  .*CEs(:) )/sum( SS_stats.mu_dist_Emp(:) );
CE_Unemp = sum( SS_stats.mu_dist_Unemp(:).*CEs(:) )/sum( SS_stats.mu_dist_Unemp(:) );

CE_I_B01   = sum( SS_stats.mu_dist_I_B01(:)  .*CEs(:) )/sum( SS_stats.mu_dist_I_B01(:) );
CE_I_B1    = sum( SS_stats.mu_dist_I_B1(:)   .*CEs(:) )/sum( SS_stats.mu_dist_I_B1(:) );
CE_I_B10   = sum( SS_stats.mu_dist_I_B10(:)  .*CEs(:) )/sum( SS_stats.mu_dist_I_B10(:) );
CE_I_Q1    = sum( SS_stats.mu_dist_I_Q1(:)   .*CEs(:) )/sum( SS_stats.mu_dist_I_Q1(:) );
CE_I_Q2    = sum( SS_stats.mu_dist_I_Q2(:)   .*CEs(:) )/sum( SS_stats.mu_dist_I_Q2(:) );
CE_I_Q3    = sum( SS_stats.mu_dist_I_Q3(:)   .*CEs(:) )/sum( SS_stats.mu_dist_I_Q3(:) );
CE_I_Q4    = sum( SS_stats.mu_dist_I_Q4(:)   .*CEs(:) )/sum( SS_stats.mu_dist_I_Q4(:) );
CE_I_Q5    = sum( SS_stats.mu_dist_I_Q5(:)   .*CEs(:) )/sum( SS_stats.mu_dist_I_Q5(:) );
CE_I_P90   = sum( SS_stats.mu_dist_I_P90(:)  .*CEs(:) )/sum( SS_stats.mu_dist_I_P90(:) );
CE_I_P99   = sum( SS_stats.mu_dist_I_P99(:)  .*CEs(:) )/sum( SS_stats.mu_dist_I_P99(:) );
CE_I_P999  = sum( SS_stats.mu_dist_I_P999(:) .*CEs(:) )/sum( SS_stats.mu_dist_I_P999(:) );

CE_Ent_Q1     = sum( (SS_stats.mu_dist_Q1(:)>0)  .*SS_stats.mu_dist_Ent(:)  .*CEs(:) )/sum( (SS_stats.mu_dist_Q1(:)>0)  .*SS_stats.mu_dist_Ent(:) );
CE_Ent_Q2     = sum( (SS_stats.mu_dist_Q2(:)>0)  .*SS_stats.mu_dist_Ent(:)  .*CEs(:) )/sum( (SS_stats.mu_dist_Q2(:)>0)  .*SS_stats.mu_dist_Ent(:) );
CE_Ent_Q3     = sum( (SS_stats.mu_dist_Q3(:)>0)  .*SS_stats.mu_dist_Ent(:)  .*CEs(:) )/sum( (SS_stats.mu_dist_Q3(:)>0)  .*SS_stats.mu_dist_Ent(:) );
CE_Ent_Q4     = sum( (SS_stats.mu_dist_Q4(:)>0)  .*SS_stats.mu_dist_Ent(:)  .*CEs(:) )/sum( (SS_stats.mu_dist_Q4(:)>0)  .*SS_stats.mu_dist_Ent(:) );
CE_Ent_Q5     = sum( (SS_stats.mu_dist_Q5(:)>0)  .*SS_stats.mu_dist_Ent(:)  .*CEs(:) )/sum( (SS_stats.mu_dist_Q5(:)>0)  .*SS_stats.mu_dist_Ent(:) );
CE_Ent_P90    = sum( (SS_stats.mu_dist_P90(:)>0) .*SS_stats.mu_dist_Ent(:)  .*CEs(:) )/sum( (SS_stats.mu_dist_P90(:)>0) .*SS_stats.mu_dist_Ent(:) );
CE_Ent_P99    = sum( (SS_stats.mu_dist_P99(:)>0) .*SS_stats.mu_dist_Ent(:)  .*CEs(:) )/sum( (SS_stats.mu_dist_P99(:)>0) .*SS_stats.mu_dist_Ent(:) );
CE_Ent_P999   = sum( (SS_stats.mu_dist_P999(:)>0).*SS_stats.mu_dist_Ent(:)  .*CEs(:) )/sum( (SS_stats.mu_dist_P999(:)>0).*SS_stats.mu_dist_Ent(:) );

CE_Emp_Q1     = sum( (SS_stats.mu_dist_Q1(:)>0)  .*SS_stats.mu_dist_Emp(:)  .*CEs(:) )/sum( (SS_stats.mu_dist_Q1(:)>0)  .*SS_stats.mu_dist_Emp(:) );
CE_Emp_Q2     = sum( (SS_stats.mu_dist_Q2(:)>0)  .*SS_stats.mu_dist_Emp(:)  .*CEs(:) )/sum( (SS_stats.mu_dist_Q2(:)>0)  .*SS_stats.mu_dist_Emp(:) );
CE_Emp_Q3     = sum( (SS_stats.mu_dist_Q3(:)>0)  .*SS_stats.mu_dist_Emp(:)  .*CEs(:) )/sum( (SS_stats.mu_dist_Q3(:)>0)  .*SS_stats.mu_dist_Emp(:) );
CE_Emp_Q4     = sum( (SS_stats.mu_dist_Q4(:)>0)  .*SS_stats.mu_dist_Emp(:)  .*CEs(:) )/sum( (SS_stats.mu_dist_Q4(:)>0)  .*SS_stats.mu_dist_Emp(:) );
CE_Emp_Q5     = sum( (SS_stats.mu_dist_Q5(:)>0)  .*SS_stats.mu_dist_Emp(:)  .*CEs(:) )/sum( (SS_stats.mu_dist_Q5(:)>0)  .*SS_stats.mu_dist_Emp(:) );
CE_Emp_P90    = sum( (SS_stats.mu_dist_P90(:)>0) .*SS_stats.mu_dist_Emp(:)  .*CEs(:) )/sum( (SS_stats.mu_dist_P90(:)>0) .*SS_stats.mu_dist_Emp(:) );
CE_Emp_P99    = sum( (SS_stats.mu_dist_P99(:)>0) .*SS_stats.mu_dist_Emp(:)  .*CEs(:) )/sum( (SS_stats.mu_dist_P99(:)>0) .*SS_stats.mu_dist_Emp(:) );
CE_Emp_P999   = sum( (SS_stats.mu_dist_P999(:)>0).*SS_stats.mu_dist_Emp(:)  .*CEs(:) )/sum( (SS_stats.mu_dist_P999(:)>0).*SS_stats.mu_dist_Emp(:) );

CE_Unemp_Q1   = sum( (SS_stats.mu_dist_Q1(:)>0)  .*SS_stats.mu_dist_Unemp(:)  .*CEs(:) )/sum( (SS_stats.mu_dist_Q1(:)>0)  .*SS_stats.mu_dist_Unemp(:) );
CE_Unemp_Q2   = sum( (SS_stats.mu_dist_Q2(:)>0)  .*SS_stats.mu_dist_Unemp(:)  .*CEs(:) )/sum( (SS_stats.mu_dist_Q2(:)>0)  .*SS_stats.mu_dist_Unemp(:) );
CE_Unemp_Q3   = sum( (SS_stats.mu_dist_Q3(:)>0)  .*SS_stats.mu_dist_Unemp(:)  .*CEs(:) )/sum( (SS_stats.mu_dist_Q3(:)>0)  .*SS_stats.mu_dist_Unemp(:) );
CE_Unemp_Q4   = sum( (SS_stats.mu_dist_Q4(:)>0)  .*SS_stats.mu_dist_Unemp(:)  .*CEs(:) )/sum( (SS_stats.mu_dist_Q4(:)>0)  .*SS_stats.mu_dist_Unemp(:) );
CE_Unemp_Q5   = sum( (SS_stats.mu_dist_Q5(:)>0)  .*SS_stats.mu_dist_Unemp(:)  .*CEs(:) )/sum( (SS_stats.mu_dist_Q5(:)>0)  .*SS_stats.mu_dist_Unemp(:) );
CE_Unemp_P90  = sum( (SS_stats.mu_dist_P90(:)>0) .*SS_stats.mu_dist_Unemp(:)  .*CEs(:) )/sum( (SS_stats.mu_dist_P90(:)>0) .*SS_stats.mu_dist_Unemp(:) );
CE_Unemp_P99  = sum( (SS_stats.mu_dist_P99(:)>0) .*SS_stats.mu_dist_Unemp(:)  .*CEs(:) )/sum( (SS_stats.mu_dist_P99(:)>0) .*SS_stats.mu_dist_Unemp(:) );
CE_Unemp_P999 = sum( (SS_stats.mu_dist_P999(:)>0).*SS_stats.mu_dist_Unemp(:)  .*CEs(:) )/sum( (SS_stats.mu_dist_P999(:)>0).*SS_stats.mu_dist_Unemp(:) );


CE_PROFIT       = sum( SS_stats.mu_dist(:)      .*CEs_PROFIT(:) )/sum( SS_stats.mu_dist(:) );
CE_PROFIT_B01   = sum( SS_stats.mu_dist_B01(:)  .*CEs_PROFIT(:) )/sum( SS_stats.mu_dist_B01(:) );
CE_PROFIT_B1    = sum( SS_stats.mu_dist_B1(:)   .*CEs_PROFIT(:) )/sum( SS_stats.mu_dist_B1(:) );
CE_PROFIT_B10   = sum( SS_stats.mu_dist_B10(:)  .*CEs_PROFIT(:) )/sum( SS_stats.mu_dist_B10(:) );
CE_PROFIT_Q1    = sum( SS_stats.mu_dist_Q1(:)   .*CEs_PROFIT(:) )/sum( SS_stats.mu_dist_Q1(:) );
CE_PROFIT_Q2    = sum( SS_stats.mu_dist_Q2(:)   .*CEs_PROFIT(:) )/sum( SS_stats.mu_dist_Q2(:) );
CE_PROFIT_Q3    = sum( SS_stats.mu_dist_Q3(:)   .*CEs_PROFIT(:) )/sum( SS_stats.mu_dist_Q3(:) );
CE_PROFIT_Q4    = sum( SS_stats.mu_dist_Q4(:)   .*CEs_PROFIT(:) )/sum( SS_stats.mu_dist_Q4(:) );
CE_PROFIT_Q5    = sum( SS_stats.mu_dist_Q5(:)   .*CEs_PROFIT(:) )/sum( SS_stats.mu_dist_Q5(:) );
CE_PROFIT_P90   = sum( SS_stats.mu_dist_P90(:)  .*CEs_PROFIT(:) )/sum( SS_stats.mu_dist_P90(:) );
CE_PROFIT_P99   = sum( SS_stats.mu_dist_P99(:)  .*CEs_PROFIT(:) )/sum( SS_stats.mu_dist_P99(:) );
CE_PROFIT_P999  = sum( SS_stats.mu_dist_P999(:) .*CEs_PROFIT(:) )/sum( SS_stats.mu_dist_P999(:) );
CE_PROFIT_Ent   = sum( SS_stats.mu_dist_Ent(:)  .*CEs_PROFIT(:) )/sum( SS_stats.mu_dist_Ent(:) );
CE_PROFIT_Emp   = sum( SS_stats.mu_dist_Emp(:)  .*CEs_PROFIT(:) )/sum( SS_stats.mu_dist_Emp(:) );
CE_PROFIT_Unemp = sum( SS_stats.mu_dist_Unemp(:).*CEs_PROFIT(:) )/sum( SS_stats.mu_dist_Unemp(:) );

CE_I_PROFIT_B01   = sum( SS_stats.mu_dist_I_B01(:)  .*CEs_PROFIT(:) )/sum( SS_stats.mu_dist_I_B01(:) );
CE_I_PROFIT_B1    = sum( SS_stats.mu_dist_I_B1(:)   .*CEs_PROFIT(:) )/sum( SS_stats.mu_dist_I_B1(:) );
CE_I_PROFIT_B10   = sum( SS_stats.mu_dist_I_B10(:)  .*CEs_PROFIT(:) )/sum( SS_stats.mu_dist_I_B10(:) );
CE_I_PROFIT_Q1    = sum( SS_stats.mu_dist_I_Q1(:)   .*CEs_PROFIT(:) )/sum( SS_stats.mu_dist_I_Q1(:) );
CE_I_PROFIT_Q2    = sum( SS_stats.mu_dist_I_Q2(:)   .*CEs_PROFIT(:) )/sum( SS_stats.mu_dist_I_Q2(:) );
CE_I_PROFIT_Q3    = sum( SS_stats.mu_dist_I_Q3(:)   .*CEs_PROFIT(:) )/sum( SS_stats.mu_dist_I_Q3(:) );
CE_I_PROFIT_Q4    = sum( SS_stats.mu_dist_I_Q4(:)   .*CEs_PROFIT(:) )/sum( SS_stats.mu_dist_I_Q4(:) );
CE_I_PROFIT_Q5    = sum( SS_stats.mu_dist_I_Q5(:)   .*CEs_PROFIT(:) )/sum( SS_stats.mu_dist_I_Q5(:) );
CE_I_PROFIT_P90   = sum( SS_stats.mu_dist_I_P90(:)  .*CEs_PROFIT(:) )/sum( SS_stats.mu_dist_I_P90(:) );
CE_I_PROFIT_P99   = sum( SS_stats.mu_dist_I_P99(:)  .*CEs_PROFIT(:) )/sum( SS_stats.mu_dist_I_P99(:) );
CE_I_PROFIT_P999  = sum( SS_stats.mu_dist_I_P999(:) .*CEs_PROFIT(:) )/sum( SS_stats.mu_dist_I_P999(:) );


CE_PROFIT_Ent_Q1     = sum( (SS_stats.mu_dist_Q1(:)>0)  .*SS_stats.mu_dist_Ent(:)  .*CEs_PROFIT(:) )/sum( (SS_stats.mu_dist_Q1(:)>0)  .*SS_stats.mu_dist_Ent(:) );
CE_PROFIT_Ent_Q2     = sum( (SS_stats.mu_dist_Q2(:)>0)  .*SS_stats.mu_dist_Ent(:)  .*CEs_PROFIT(:) )/sum( (SS_stats.mu_dist_Q2(:)>0)  .*SS_stats.mu_dist_Ent(:) );
CE_PROFIT_Ent_Q3     = sum( (SS_stats.mu_dist_Q3(:)>0)  .*SS_stats.mu_dist_Ent(:)  .*CEs_PROFIT(:) )/sum( (SS_stats.mu_dist_Q3(:)>0)  .*SS_stats.mu_dist_Ent(:) );
CE_PROFIT_Ent_Q4     = sum( (SS_stats.mu_dist_Q4(:)>0)  .*SS_stats.mu_dist_Ent(:)  .*CEs_PROFIT(:) )/sum( (SS_stats.mu_dist_Q4(:)>0)  .*SS_stats.mu_dist_Ent(:) );
CE_PROFIT_Ent_Q5     = sum( (SS_stats.mu_dist_Q5(:)>0)  .*SS_stats.mu_dist_Ent(:)  .*CEs_PROFIT(:) )/sum( (SS_stats.mu_dist_Q5(:)>0)  .*SS_stats.mu_dist_Ent(:) );
CE_PROFIT_Ent_P90    = sum( (SS_stats.mu_dist_P90(:)>0) .*SS_stats.mu_dist_Ent(:)  .*CEs_PROFIT(:) )/sum( (SS_stats.mu_dist_P90(:)>0) .*SS_stats.mu_dist_Ent(:) );
CE_PROFIT_Ent_P99    = sum( (SS_stats.mu_dist_P99(:)>0) .*SS_stats.mu_dist_Ent(:)  .*CEs_PROFIT(:) )/sum( (SS_stats.mu_dist_P99(:)>0) .*SS_stats.mu_dist_Ent(:) );
CE_PROFIT_Ent_P999   = sum( (SS_stats.mu_dist_P999(:)>0).*SS_stats.mu_dist_Ent(:)  .*CEs_PROFIT(:) )/sum( (SS_stats.mu_dist_P999(:)>0).*SS_stats.mu_dist_Ent(:) );

CE_PROFIT_Emp_Q1     = sum( (SS_stats.mu_dist_Q1(:)>0)  .*SS_stats.mu_dist_Emp(:)  .*CEs_PROFIT(:) )/sum( (SS_stats.mu_dist_Q1(:)>0)  .*SS_stats.mu_dist_Emp(:) );
CE_PROFIT_Emp_Q2     = sum( (SS_stats.mu_dist_Q2(:)>0)  .*SS_stats.mu_dist_Emp(:)  .*CEs_PROFIT(:) )/sum( (SS_stats.mu_dist_Q2(:)>0)  .*SS_stats.mu_dist_Emp(:) );
CE_PROFIT_Emp_Q3     = sum( (SS_stats.mu_dist_Q3(:)>0)  .*SS_stats.mu_dist_Emp(:)  .*CEs_PROFIT(:) )/sum( (SS_stats.mu_dist_Q3(:)>0)  .*SS_stats.mu_dist_Emp(:) );
CE_PROFIT_Emp_Q4     = sum( (SS_stats.mu_dist_Q4(:)>0)  .*SS_stats.mu_dist_Emp(:)  .*CEs_PROFIT(:) )/sum( (SS_stats.mu_dist_Q4(:)>0)  .*SS_stats.mu_dist_Emp(:) );
CE_PROFIT_Emp_Q5     = sum( (SS_stats.mu_dist_Q5(:)>0)  .*SS_stats.mu_dist_Emp(:)  .*CEs_PROFIT(:) )/sum( (SS_stats.mu_dist_Q5(:)>0)  .*SS_stats.mu_dist_Emp(:) );
CE_PROFIT_Emp_P90    = sum( (SS_stats.mu_dist_P90(:)>0) .*SS_stats.mu_dist_Emp(:)  .*CEs_PROFIT(:) )/sum( (SS_stats.mu_dist_P90(:)>0) .*SS_stats.mu_dist_Emp(:) );
CE_PROFIT_Emp_P99    = sum( (SS_stats.mu_dist_P99(:)>0) .*SS_stats.mu_dist_Emp(:)  .*CEs_PROFIT(:) )/sum( (SS_stats.mu_dist_P99(:)>0) .*SS_stats.mu_dist_Emp(:) );
CE_PROFIT_Emp_P999   = sum( (SS_stats.mu_dist_P999(:)>0).*SS_stats.mu_dist_Emp(:)  .*CEs_PROFIT(:) )/sum( (SS_stats.mu_dist_P999(:)>0).*SS_stats.mu_dist_Emp(:) );

CE_PROFIT_Unemp_Q1   = sum( (SS_stats.mu_dist_Q1(:)>0)  .*SS_stats.mu_dist_Unemp(:)  .*CEs_PROFIT(:) )/sum( (SS_stats.mu_dist_Q1(:)>0)  .*SS_stats.mu_dist_Unemp(:) );
CE_PROFIT_Unemp_Q2   = sum( (SS_stats.mu_dist_Q2(:)>0)  .*SS_stats.mu_dist_Unemp(:)  .*CEs_PROFIT(:) )/sum( (SS_stats.mu_dist_Q2(:)>0)  .*SS_stats.mu_dist_Unemp(:) );
CE_PROFIT_Unemp_Q3   = sum( (SS_stats.mu_dist_Q3(:)>0)  .*SS_stats.mu_dist_Unemp(:)  .*CEs_PROFIT(:) )/sum( (SS_stats.mu_dist_Q3(:)>0)  .*SS_stats.mu_dist_Unemp(:) );
CE_PROFIT_Unemp_Q4   = sum( (SS_stats.mu_dist_Q4(:)>0)  .*SS_stats.mu_dist_Unemp(:)  .*CEs_PROFIT(:) )/sum( (SS_stats.mu_dist_Q4(:)>0)  .*SS_stats.mu_dist_Unemp(:) );
CE_PROFIT_Unemp_Q5   = sum( (SS_stats.mu_dist_Q5(:)>0)  .*SS_stats.mu_dist_Unemp(:)  .*CEs_PROFIT(:) )/sum( (SS_stats.mu_dist_Q5(:)>0)  .*SS_stats.mu_dist_Unemp(:) );
CE_PROFIT_Unemp_P90  = sum( (SS_stats.mu_dist_P90(:)>0) .*SS_stats.mu_dist_Unemp(:)  .*CEs_PROFIT(:) )/sum( (SS_stats.mu_dist_P90(:)>0) .*SS_stats.mu_dist_Unemp(:) );
CE_PROFIT_Unemp_P99  = sum( (SS_stats.mu_dist_P99(:)>0) .*SS_stats.mu_dist_Unemp(:)  .*CEs_PROFIT(:) )/sum( (SS_stats.mu_dist_P99(:)>0) .*SS_stats.mu_dist_Unemp(:) );
CE_PROFIT_Unemp_P999 = sum( (SS_stats.mu_dist_P999(:)>0).*SS_stats.mu_dist_Unemp(:)  .*CEs_PROFIT(:) )/sum( (SS_stats.mu_dist_P999(:)>0).*SS_stats.mu_dist_Unemp(:) );



CE_RR       = sum( SS_stats.mu_dist(:)      .*CEs_RR(:) )/sum( SS_stats.mu_dist(:) );
CE_RR_B01   = sum( SS_stats.mu_dist_B01(:)  .*CEs_RR(:) )/sum( SS_stats.mu_dist_B01(:) );
CE_RR_B1    = sum( SS_stats.mu_dist_B1(:)   .*CEs_RR(:) )/sum( SS_stats.mu_dist_B1(:) );
CE_RR_B10   = sum( SS_stats.mu_dist_B10(:)  .*CEs_RR(:) )/sum( SS_stats.mu_dist_B10(:) );
CE_RR_Q1    = sum( SS_stats.mu_dist_Q1(:)   .*CEs_RR(:) )/sum( SS_stats.mu_dist_Q1(:) );
CE_RR_Q2    = sum( SS_stats.mu_dist_Q2(:)   .*CEs_RR(:) )/sum( SS_stats.mu_dist_Q2(:) );
CE_RR_Q3    = sum( SS_stats.mu_dist_Q3(:)   .*CEs_RR(:) )/sum( SS_stats.mu_dist_Q3(:) );
CE_RR_Q4    = sum( SS_stats.mu_dist_Q4(:)   .*CEs_RR(:) )/sum( SS_stats.mu_dist_Q4(:) );
CE_RR_Q5    = sum( SS_stats.mu_dist_Q5(:)   .*CEs_RR(:) )/sum( SS_stats.mu_dist_Q5(:) );
CE_RR_P90   = sum( SS_stats.mu_dist_P90(:)  .*CEs_RR(:) )/sum( SS_stats.mu_dist_P90(:) );
CE_RR_P99   = sum( SS_stats.mu_dist_P99(:)  .*CEs_RR(:) )/sum( SS_stats.mu_dist_P99(:) );
CE_RR_P999  = sum( SS_stats.mu_dist_P999(:) .*CEs_RR(:) )/sum( SS_stats.mu_dist_P999(:) );
CE_RR_Ent   = sum( SS_stats.mu_dist_Ent(:)  .*CEs_RR(:) )/sum( SS_stats.mu_dist_Ent(:) );
CE_RR_Emp   = sum( SS_stats.mu_dist_Emp(:)  .*CEs_RR(:) )/sum( SS_stats.mu_dist_Emp(:) );
CE_RR_Unemp = sum( SS_stats.mu_dist_Unemp(:).*CEs_RR(:) )/sum( SS_stats.mu_dist_Unemp(:) );

CE_I_RR_B01   = sum( SS_stats.mu_dist_I_B01(:)  .*CEs_RR(:) )/sum( SS_stats.mu_dist_I_B01(:) );
CE_I_RR_B1    = sum( SS_stats.mu_dist_I_B1(:)   .*CEs_RR(:) )/sum( SS_stats.mu_dist_I_B1(:) );
CE_I_RR_B10   = sum( SS_stats.mu_dist_I_B10(:)  .*CEs_RR(:) )/sum( SS_stats.mu_dist_I_B10(:) );
CE_I_RR_Q1    = sum( SS_stats.mu_dist_I_Q1(:)   .*CEs_RR(:) )/sum( SS_stats.mu_dist_I_Q1(:) );
CE_I_RR_Q2    = sum( SS_stats.mu_dist_I_Q2(:)   .*CEs_RR(:) )/sum( SS_stats.mu_dist_I_Q2(:) );
CE_I_RR_Q3    = sum( SS_stats.mu_dist_I_Q3(:)   .*CEs_RR(:) )/sum( SS_stats.mu_dist_I_Q3(:) );
CE_I_RR_Q4    = sum( SS_stats.mu_dist_I_Q4(:)   .*CEs_RR(:) )/sum( SS_stats.mu_dist_I_Q4(:) );
CE_I_RR_Q5    = sum( SS_stats.mu_dist_I_Q5(:)   .*CEs_RR(:) )/sum( SS_stats.mu_dist_I_Q5(:) );
CE_I_RR_P90   = sum( SS_stats.mu_dist_I_P90(:)  .*CEs_RR(:) )/sum( SS_stats.mu_dist_I_P90(:) );
CE_I_RR_P99   = sum( SS_stats.mu_dist_I_P99(:)  .*CEs_RR(:) )/sum( SS_stats.mu_dist_I_P99(:) );
CE_I_RR_P999  = sum( SS_stats.mu_dist_I_P999(:) .*CEs_RR(:) )/sum( SS_stats.mu_dist_I_P999(:) );


CE_RR_Ent_Q1     = sum( (SS_stats.mu_dist_Q1(:)>0)  .*SS_stats.mu_dist_Ent(:)  .*CEs_RR(:) )/sum( (SS_stats.mu_dist_Q1(:)>0)  .*SS_stats.mu_dist_Ent(:) );
CE_RR_Ent_Q2     = sum( (SS_stats.mu_dist_Q2(:)>0)  .*SS_stats.mu_dist_Ent(:)  .*CEs_RR(:) )/sum( (SS_stats.mu_dist_Q2(:)>0)  .*SS_stats.mu_dist_Ent(:) );
CE_RR_Ent_Q3     = sum( (SS_stats.mu_dist_Q3(:)>0)  .*SS_stats.mu_dist_Ent(:)  .*CEs_RR(:) )/sum( (SS_stats.mu_dist_Q3(:)>0)  .*SS_stats.mu_dist_Ent(:) );
CE_RR_Ent_Q4     = sum( (SS_stats.mu_dist_Q4(:)>0)  .*SS_stats.mu_dist_Ent(:)  .*CEs_RR(:) )/sum( (SS_stats.mu_dist_Q4(:)>0)  .*SS_stats.mu_dist_Ent(:) );
CE_RR_Ent_Q5     = sum( (SS_stats.mu_dist_Q5(:)>0)  .*SS_stats.mu_dist_Ent(:)  .*CEs_RR(:) )/sum( (SS_stats.mu_dist_Q5(:)>0)  .*SS_stats.mu_dist_Ent(:) );
CE_RR_Ent_P90    = sum( (SS_stats.mu_dist_P90(:)>0) .*SS_stats.mu_dist_Ent(:)  .*CEs_RR(:) )/sum( (SS_stats.mu_dist_P90(:)>0) .*SS_stats.mu_dist_Ent(:) );
CE_RR_Ent_P99    = sum( (SS_stats.mu_dist_P99(:)>0) .*SS_stats.mu_dist_Ent(:)  .*CEs_RR(:) )/sum( (SS_stats.mu_dist_P99(:)>0) .*SS_stats.mu_dist_Ent(:) );
CE_RR_Ent_P999   = sum( (SS_stats.mu_dist_P999(:)>0).*SS_stats.mu_dist_Ent(:)  .*CEs_RR(:) )/sum( (SS_stats.mu_dist_P999(:)>0).*SS_stats.mu_dist_Ent(:) );

CE_RR_Emp_Q1     = sum( (SS_stats.mu_dist_Q1(:)>0)  .*SS_stats.mu_dist_Emp(:)  .*CEs_RR(:) )/sum( (SS_stats.mu_dist_Q1(:)>0)  .*SS_stats.mu_dist_Emp(:) );
CE_RR_Emp_Q2     = sum( (SS_stats.mu_dist_Q2(:)>0)  .*SS_stats.mu_dist_Emp(:)  .*CEs_RR(:) )/sum( (SS_stats.mu_dist_Q2(:)>0)  .*SS_stats.mu_dist_Emp(:) );
CE_RR_Emp_Q3     = sum( (SS_stats.mu_dist_Q3(:)>0)  .*SS_stats.mu_dist_Emp(:)  .*CEs_RR(:) )/sum( (SS_stats.mu_dist_Q3(:)>0)  .*SS_stats.mu_dist_Emp(:) );
CE_RR_Emp_Q4     = sum( (SS_stats.mu_dist_Q4(:)>0)  .*SS_stats.mu_dist_Emp(:)  .*CEs_RR(:) )/sum( (SS_stats.mu_dist_Q4(:)>0)  .*SS_stats.mu_dist_Emp(:) );
CE_RR_Emp_Q5     = sum( (SS_stats.mu_dist_Q5(:)>0)  .*SS_stats.mu_dist_Emp(:)  .*CEs_RR(:) )/sum( (SS_stats.mu_dist_Q5(:)>0)  .*SS_stats.mu_dist_Emp(:) );
CE_RR_Emp_P90    = sum( (SS_stats.mu_dist_P90(:)>0) .*SS_stats.mu_dist_Emp(:)  .*CEs_RR(:) )/sum( (SS_stats.mu_dist_P90(:)>0) .*SS_stats.mu_dist_Emp(:) );
CE_RR_Emp_P99    = sum( (SS_stats.mu_dist_P99(:)>0) .*SS_stats.mu_dist_Emp(:)  .*CEs_RR(:) )/sum( (SS_stats.mu_dist_P99(:)>0) .*SS_stats.mu_dist_Emp(:) );
CE_RR_Emp_P999   = sum( (SS_stats.mu_dist_P999(:)>0).*SS_stats.mu_dist_Emp(:)  .*CEs_RR(:) )/sum( (SS_stats.mu_dist_P999(:)>0).*SS_stats.mu_dist_Emp(:) );

CE_RR_Unemp_Q1   = sum( (SS_stats.mu_dist_Q1(:)>0)  .*SS_stats.mu_dist_Unemp(:)  .*CEs_RR(:) )/sum( (SS_stats.mu_dist_Q1(:)>0)  .*SS_stats.mu_dist_Unemp(:) );
CE_RR_Unemp_Q2   = sum( (SS_stats.mu_dist_Q2(:)>0)  .*SS_stats.mu_dist_Unemp(:)  .*CEs_RR(:) )/sum( (SS_stats.mu_dist_Q2(:)>0)  .*SS_stats.mu_dist_Unemp(:) );
CE_RR_Unemp_Q3   = sum( (SS_stats.mu_dist_Q3(:)>0)  .*SS_stats.mu_dist_Unemp(:)  .*CEs_RR(:) )/sum( (SS_stats.mu_dist_Q3(:)>0)  .*SS_stats.mu_dist_Unemp(:) );
CE_RR_Unemp_Q4   = sum( (SS_stats.mu_dist_Q4(:)>0)  .*SS_stats.mu_dist_Unemp(:)  .*CEs_RR(:) )/sum( (SS_stats.mu_dist_Q4(:)>0)  .*SS_stats.mu_dist_Unemp(:) );
CE_RR_Unemp_Q5   = sum( (SS_stats.mu_dist_Q5(:)>0)  .*SS_stats.mu_dist_Unemp(:)  .*CEs_RR(:) )/sum( (SS_stats.mu_dist_Q5(:)>0)  .*SS_stats.mu_dist_Unemp(:) );
CE_RR_Unemp_P90  = sum( (SS_stats.mu_dist_P90(:)>0) .*SS_stats.mu_dist_Unemp(:)  .*CEs_RR(:) )/sum( (SS_stats.mu_dist_P90(:)>0) .*SS_stats.mu_dist_Unemp(:) );
CE_RR_Unemp_P99  = sum( (SS_stats.mu_dist_P99(:)>0) .*SS_stats.mu_dist_Unemp(:)  .*CEs_RR(:) )/sum( (SS_stats.mu_dist_P99(:)>0) .*SS_stats.mu_dist_Unemp(:) );
CE_RR_Unemp_P999 = sum( (SS_stats.mu_dist_P999(:)>0).*SS_stats.mu_dist_Unemp(:)  .*CEs_RR(:) )/sum( (SS_stats.mu_dist_P999(:)>0).*SS_stats.mu_dist_Unemp(:) );

CE_W       = sum( SS_stats.mu_dist(:)      .*CEs_W(:) )/sum( SS_stats.mu_dist(:) );
CE_W_B01   = sum( SS_stats.mu_dist_B01(:)  .*CEs_W(:) )/sum( SS_stats.mu_dist_B01(:) );
CE_W_B1    = sum( SS_stats.mu_dist_B1(:)   .*CEs_W(:) )/sum( SS_stats.mu_dist_B1(:) );
CE_W_B10   = sum( SS_stats.mu_dist_B10(:)  .*CEs_W(:) )/sum( SS_stats.mu_dist_B10(:) );
CE_W_Q1    = sum( SS_stats.mu_dist_Q1(:)   .*CEs_W(:) )/sum( SS_stats.mu_dist_Q1(:) );
CE_W_Q2    = sum( SS_stats.mu_dist_Q2(:)   .*CEs_W(:) )/sum( SS_stats.mu_dist_Q2(:) );
CE_W_Q3    = sum( SS_stats.mu_dist_Q3(:)   .*CEs_W(:) )/sum( SS_stats.mu_dist_Q3(:) );
CE_W_Q4    = sum( SS_stats.mu_dist_Q4(:)   .*CEs_W(:) )/sum( SS_stats.mu_dist_Q4(:) );
CE_W_Q5    = sum( SS_stats.mu_dist_Q5(:)   .*CEs_W(:) )/sum( SS_stats.mu_dist_Q5(:) );
CE_W_P90   = sum( SS_stats.mu_dist_P90(:)  .*CEs_W(:) )/sum( SS_stats.mu_dist_P90(:) );
CE_W_P99   = sum( SS_stats.mu_dist_P99(:)  .*CEs_W(:) )/sum( SS_stats.mu_dist_P99(:) );
CE_W_P999  = sum( SS_stats.mu_dist_P999(:) .*CEs_W(:) )/sum( SS_stats.mu_dist_P999(:) );
CE_W_Ent   = sum( SS_stats.mu_dist_Ent(:)  .*CEs_W(:) )/sum( SS_stats.mu_dist_Ent(:) );
CE_W_Emp   = sum( SS_stats.mu_dist_Emp(:)  .*CEs_W(:) )/sum( SS_stats.mu_dist_Emp(:) );
CE_W_Unemp = sum( SS_stats.mu_dist_Unemp(:).*CEs_W(:) )/sum( SS_stats.mu_dist_Unemp(:) );

CE_I_W_B01   = sum( SS_stats.mu_dist_I_B01(:)  .*CEs_W(:) )/sum( SS_stats.mu_dist_I_B01(:) );
CE_I_W_B1    = sum( SS_stats.mu_dist_I_B1(:)   .*CEs_W(:) )/sum( SS_stats.mu_dist_I_B1(:) );
CE_I_W_B10   = sum( SS_stats.mu_dist_I_B10(:)  .*CEs_W(:) )/sum( SS_stats.mu_dist_I_B10(:) );
CE_I_W_Q1    = sum( SS_stats.mu_dist_I_Q1(:)   .*CEs_W(:) )/sum( SS_stats.mu_dist_I_Q1(:) );
CE_I_W_Q2    = sum( SS_stats.mu_dist_I_Q2(:)   .*CEs_W(:) )/sum( SS_stats.mu_dist_I_Q2(:) );
CE_I_W_Q3    = sum( SS_stats.mu_dist_I_Q3(:)   .*CEs_W(:) )/sum( SS_stats.mu_dist_I_Q3(:) );
CE_I_W_Q4    = sum( SS_stats.mu_dist_I_Q4(:)   .*CEs_W(:) )/sum( SS_stats.mu_dist_I_Q4(:) );
CE_I_W_Q5    = sum( SS_stats.mu_dist_I_Q5(:)   .*CEs_W(:) )/sum( SS_stats.mu_dist_I_Q5(:) );
CE_I_W_P90   = sum( SS_stats.mu_dist_I_P90(:)  .*CEs_W(:) )/sum( SS_stats.mu_dist_I_P90(:) );
CE_I_W_P99   = sum( SS_stats.mu_dist_I_P99(:)  .*CEs_W(:) )/sum( SS_stats.mu_dist_I_P99(:) );
CE_I_W_P999  = sum( SS_stats.mu_dist_I_P999(:) .*CEs_W(:) )/sum( SS_stats.mu_dist_I_P999(:) );


CE_W_Ent_Q1     = sum( (SS_stats.mu_dist_Q1(:)>0)  .*SS_stats.mu_dist_Ent(:)  .*CEs_W(:) )/sum( (SS_stats.mu_dist_Q1(:)>0)  .*SS_stats.mu_dist_Ent(:) );
CE_W_Ent_Q2     = sum( (SS_stats.mu_dist_Q2(:)>0)  .*SS_stats.mu_dist_Ent(:)  .*CEs_W(:) )/sum( (SS_stats.mu_dist_Q2(:)>0)  .*SS_stats.mu_dist_Ent(:) );
CE_W_Ent_Q3     = sum( (SS_stats.mu_dist_Q3(:)>0)  .*SS_stats.mu_dist_Ent(:)  .*CEs_W(:) )/sum( (SS_stats.mu_dist_Q3(:)>0)  .*SS_stats.mu_dist_Ent(:) );
CE_W_Ent_Q4     = sum( (SS_stats.mu_dist_Q4(:)>0)  .*SS_stats.mu_dist_Ent(:)  .*CEs_W(:) )/sum( (SS_stats.mu_dist_Q4(:)>0)  .*SS_stats.mu_dist_Ent(:) );
CE_W_Ent_Q5     = sum( (SS_stats.mu_dist_Q5(:)>0)  .*SS_stats.mu_dist_Ent(:)  .*CEs_W(:) )/sum( (SS_stats.mu_dist_Q5(:)>0)  .*SS_stats.mu_dist_Ent(:) );
CE_W_Ent_P90    = sum( (SS_stats.mu_dist_P90(:)>0) .*SS_stats.mu_dist_Ent(:)  .*CEs_W(:) )/sum( (SS_stats.mu_dist_P90(:)>0) .*SS_stats.mu_dist_Ent(:) );
CE_W_Ent_P99    = sum( (SS_stats.mu_dist_P99(:)>0) .*SS_stats.mu_dist_Ent(:)  .*CEs_W(:) )/sum( (SS_stats.mu_dist_P99(:)>0) .*SS_stats.mu_dist_Ent(:) );
CE_W_Ent_P999   = sum( (SS_stats.mu_dist_P999(:)>0).*SS_stats.mu_dist_Ent(:)  .*CEs_W(:) )/sum( (SS_stats.mu_dist_P999(:)>0).*SS_stats.mu_dist_Ent(:) );

CE_W_Emp_Q1     = sum( (SS_stats.mu_dist_Q1(:)>0)  .*SS_stats.mu_dist_Emp(:)  .*CEs_W(:) )/sum( (SS_stats.mu_dist_Q1(:)>0)  .*SS_stats.mu_dist_Emp(:) );
CE_W_Emp_Q2     = sum( (SS_stats.mu_dist_Q2(:)>0)  .*SS_stats.mu_dist_Emp(:)  .*CEs_W(:) )/sum( (SS_stats.mu_dist_Q2(:)>0)  .*SS_stats.mu_dist_Emp(:) );
CE_W_Emp_Q3     = sum( (SS_stats.mu_dist_Q3(:)>0)  .*SS_stats.mu_dist_Emp(:)  .*CEs_W(:) )/sum( (SS_stats.mu_dist_Q3(:)>0)  .*SS_stats.mu_dist_Emp(:) );
CE_W_Emp_Q4     = sum( (SS_stats.mu_dist_Q4(:)>0)  .*SS_stats.mu_dist_Emp(:)  .*CEs_W(:) )/sum( (SS_stats.mu_dist_Q4(:)>0)  .*SS_stats.mu_dist_Emp(:) );
CE_W_Emp_Q5     = sum( (SS_stats.mu_dist_Q5(:)>0)  .*SS_stats.mu_dist_Emp(:)  .*CEs_W(:) )/sum( (SS_stats.mu_dist_Q5(:)>0)  .*SS_stats.mu_dist_Emp(:) );
CE_W_Emp_P90    = sum( (SS_stats.mu_dist_P90(:)>0) .*SS_stats.mu_dist_Emp(:)  .*CEs_W(:) )/sum( (SS_stats.mu_dist_P90(:)>0) .*SS_stats.mu_dist_Emp(:) );
CE_W_Emp_P99    = sum( (SS_stats.mu_dist_P99(:)>0) .*SS_stats.mu_dist_Emp(:)  .*CEs_W(:) )/sum( (SS_stats.mu_dist_P99(:)>0) .*SS_stats.mu_dist_Emp(:) );
CE_W_Emp_P999   = sum( (SS_stats.mu_dist_P999(:)>0).*SS_stats.mu_dist_Emp(:)  .*CEs_W(:) )/sum( (SS_stats.mu_dist_P999(:)>0).*SS_stats.mu_dist_Emp(:) );

CE_W_Unemp_Q1   = sum( (SS_stats.mu_dist_Q1(:)>0)  .*SS_stats.mu_dist_Unemp(:)  .*CEs_W(:) )/sum( (SS_stats.mu_dist_Q1(:)>0)  .*SS_stats.mu_dist_Unemp(:) );
CE_W_Unemp_Q2   = sum( (SS_stats.mu_dist_Q2(:)>0)  .*SS_stats.mu_dist_Unemp(:)  .*CEs_W(:) )/sum( (SS_stats.mu_dist_Q2(:)>0)  .*SS_stats.mu_dist_Unemp(:) );
CE_W_Unemp_Q3   = sum( (SS_stats.mu_dist_Q3(:)>0)  .*SS_stats.mu_dist_Unemp(:)  .*CEs_W(:) )/sum( (SS_stats.mu_dist_Q3(:)>0)  .*SS_stats.mu_dist_Unemp(:) );
CE_W_Unemp_Q4   = sum( (SS_stats.mu_dist_Q4(:)>0)  .*SS_stats.mu_dist_Unemp(:)  .*CEs_W(:) )/sum( (SS_stats.mu_dist_Q4(:)>0)  .*SS_stats.mu_dist_Unemp(:) );
CE_W_Unemp_Q5   = sum( (SS_stats.mu_dist_Q5(:)>0)  .*SS_stats.mu_dist_Unemp(:)  .*CEs_W(:) )/sum( (SS_stats.mu_dist_Q5(:)>0)  .*SS_stats.mu_dist_Unemp(:) );
CE_W_Unemp_P90  = sum( (SS_stats.mu_dist_P90(:)>0) .*SS_stats.mu_dist_Unemp(:)  .*CEs_W(:) )/sum( (SS_stats.mu_dist_P90(:)>0) .*SS_stats.mu_dist_Unemp(:) );
CE_W_Unemp_P99  = sum( (SS_stats.mu_dist_P99(:)>0) .*SS_stats.mu_dist_Unemp(:)  .*CEs_W(:) )/sum( (SS_stats.mu_dist_P99(:)>0) .*SS_stats.mu_dist_Unemp(:) );
CE_W_Unemp_P999 = sum( (SS_stats.mu_dist_P999(:)>0).*SS_stats.mu_dist_Unemp(:)  .*CEs_W(:) )/sum( (SS_stats.mu_dist_P999(:)>0).*SS_stats.mu_dist_Unemp(:) );


CE_f       = sum( SS_stats.mu_dist(:)      .*CEs_f(:) )/sum( SS_stats.mu_dist(:) );
CE_f_B01   = sum( SS_stats.mu_dist_B01(:)  .*CEs_f(:) )/sum( SS_stats.mu_dist_B01(:) );
CE_f_B1    = sum( SS_stats.mu_dist_B1(:)   .*CEs_f(:) )/sum( SS_stats.mu_dist_B1(:) );
CE_f_B10   = sum( SS_stats.mu_dist_B10(:)  .*CEs_f(:) )/sum( SS_stats.mu_dist_B10(:) );
CE_f_Q1    = sum( SS_stats.mu_dist_Q1(:)   .*CEs_f(:) )/sum( SS_stats.mu_dist_Q1(:) );
CE_f_Q2    = sum( SS_stats.mu_dist_Q2(:)   .*CEs_f(:) )/sum( SS_stats.mu_dist_Q2(:) );
CE_f_Q3    = sum( SS_stats.mu_dist_Q3(:)   .*CEs_f(:) )/sum( SS_stats.mu_dist_Q3(:) );
CE_f_Q4    = sum( SS_stats.mu_dist_Q4(:)   .*CEs_f(:) )/sum( SS_stats.mu_dist_Q4(:) );
CE_f_Q5    = sum( SS_stats.mu_dist_Q5(:)   .*CEs_f(:) )/sum( SS_stats.mu_dist_Q5(:) );
CE_f_P90   = sum( SS_stats.mu_dist_P90(:)  .*CEs_f(:) )/sum( SS_stats.mu_dist_P90(:) );
CE_f_P99   = sum( SS_stats.mu_dist_P99(:)  .*CEs_f(:) )/sum( SS_stats.mu_dist_P99(:) );
CE_f_P999  = sum( SS_stats.mu_dist_P999(:) .*CEs_f(:) )/sum( SS_stats.mu_dist_P999(:) );
CE_f_Ent   = sum( SS_stats.mu_dist_Ent(:)  .*CEs_f(:) )/sum( SS_stats.mu_dist_Ent(:) );
CE_f_Emp   = sum( SS_stats.mu_dist_Emp(:)  .*CEs_f(:) )/sum( SS_stats.mu_dist_Emp(:) );
CE_f_Unemp = sum( SS_stats.mu_dist_Unemp(:).*CEs_f(:) )/sum( SS_stats.mu_dist_Unemp(:) );

CE_I_f_B01   = sum( SS_stats.mu_dist_I_B01(:)  .*CEs_f(:) )/sum( SS_stats.mu_dist_I_B01(:) );
CE_I_f_B1    = sum( SS_stats.mu_dist_I_B1(:)   .*CEs_f(:) )/sum( SS_stats.mu_dist_I_B1(:) );
CE_I_f_B10   = sum( SS_stats.mu_dist_I_B10(:)  .*CEs_f(:) )/sum( SS_stats.mu_dist_I_B10(:) );
CE_I_f_Q1    = sum( SS_stats.mu_dist_I_Q1(:)   .*CEs_f(:) )/sum( SS_stats.mu_dist_I_Q1(:) );
CE_I_f_Q2    = sum( SS_stats.mu_dist_I_Q2(:)   .*CEs_f(:) )/sum( SS_stats.mu_dist_I_Q2(:) );
CE_I_f_Q3    = sum( SS_stats.mu_dist_I_Q3(:)   .*CEs_f(:) )/sum( SS_stats.mu_dist_I_Q3(:) );
CE_I_f_Q4    = sum( SS_stats.mu_dist_I_Q4(:)   .*CEs_f(:) )/sum( SS_stats.mu_dist_I_Q4(:) );
CE_I_f_Q5    = sum( SS_stats.mu_dist_I_Q5(:)   .*CEs_f(:) )/sum( SS_stats.mu_dist_I_Q5(:) );
CE_I_f_P90   = sum( SS_stats.mu_dist_I_P90(:)  .*CEs_f(:) )/sum( SS_stats.mu_dist_I_P90(:) );
CE_I_f_P99   = sum( SS_stats.mu_dist_I_P99(:)  .*CEs_f(:) )/sum( SS_stats.mu_dist_I_P99(:) );
CE_I_f_P999  = sum( SS_stats.mu_dist_I_P999(:) .*CEs_f(:) )/sum( SS_stats.mu_dist_I_P999(:) );

CE_f_Ent_Q1     = sum( (SS_stats.mu_dist_Q1(:)>0)  .*SS_stats.mu_dist_Ent(:)  .*CEs_f(:) )/sum( (SS_stats.mu_dist_Q1(:)>0)  .*SS_stats.mu_dist_Ent(:) );
CE_f_Ent_Q2     = sum( (SS_stats.mu_dist_Q2(:)>0)  .*SS_stats.mu_dist_Ent(:)  .*CEs_f(:) )/sum( (SS_stats.mu_dist_Q2(:)>0)  .*SS_stats.mu_dist_Ent(:) );
CE_f_Ent_Q3     = sum( (SS_stats.mu_dist_Q3(:)>0)  .*SS_stats.mu_dist_Ent(:)  .*CEs_f(:) )/sum( (SS_stats.mu_dist_Q3(:)>0)  .*SS_stats.mu_dist_Ent(:) );
CE_f_Ent_Q4     = sum( (SS_stats.mu_dist_Q4(:)>0)  .*SS_stats.mu_dist_Ent(:)  .*CEs_f(:) )/sum( (SS_stats.mu_dist_Q4(:)>0)  .*SS_stats.mu_dist_Ent(:) );
CE_f_Ent_Q5     = sum( (SS_stats.mu_dist_Q5(:)>0)  .*SS_stats.mu_dist_Ent(:)  .*CEs_f(:) )/sum( (SS_stats.mu_dist_Q5(:)>0)  .*SS_stats.mu_dist_Ent(:) );
CE_f_Ent_P90    = sum( (SS_stats.mu_dist_P90(:)>0) .*SS_stats.mu_dist_Ent(:)  .*CEs_f(:) )/sum( (SS_stats.mu_dist_P90(:)>0) .*SS_stats.mu_dist_Ent(:) );
CE_f_Ent_P99    = sum( (SS_stats.mu_dist_P99(:)>0) .*SS_stats.mu_dist_Ent(:)  .*CEs_f(:) )/sum( (SS_stats.mu_dist_P99(:)>0) .*SS_stats.mu_dist_Ent(:) );
CE_f_Ent_P999   = sum( (SS_stats.mu_dist_P999(:)>0).*SS_stats.mu_dist_Ent(:)  .*CEs_f(:) )/sum( (SS_stats.mu_dist_P999(:)>0).*SS_stats.mu_dist_Ent(:) );

CE_f_Emp_Q1     = sum( (SS_stats.mu_dist_Q1(:)>0)  .*SS_stats.mu_dist_Emp(:)  .*CEs_f(:) )/sum( (SS_stats.mu_dist_Q1(:)>0)  .*SS_stats.mu_dist_Emp(:) );
CE_f_Emp_Q2     = sum( (SS_stats.mu_dist_Q2(:)>0)  .*SS_stats.mu_dist_Emp(:)  .*CEs_f(:) )/sum( (SS_stats.mu_dist_Q2(:)>0)  .*SS_stats.mu_dist_Emp(:) );
CE_f_Emp_Q3     = sum( (SS_stats.mu_dist_Q3(:)>0)  .*SS_stats.mu_dist_Emp(:)  .*CEs_f(:) )/sum( (SS_stats.mu_dist_Q3(:)>0)  .*SS_stats.mu_dist_Emp(:) );
CE_f_Emp_Q4     = sum( (SS_stats.mu_dist_Q4(:)>0)  .*SS_stats.mu_dist_Emp(:)  .*CEs_f(:) )/sum( (SS_stats.mu_dist_Q4(:)>0)  .*SS_stats.mu_dist_Emp(:) );
CE_f_Emp_Q5     = sum( (SS_stats.mu_dist_Q5(:)>0)  .*SS_stats.mu_dist_Emp(:)  .*CEs_f(:) )/sum( (SS_stats.mu_dist_Q5(:)>0)  .*SS_stats.mu_dist_Emp(:) );
CE_f_Emp_P90    = sum( (SS_stats.mu_dist_P90(:)>0) .*SS_stats.mu_dist_Emp(:)  .*CEs_f(:) )/sum( (SS_stats.mu_dist_P90(:)>0) .*SS_stats.mu_dist_Emp(:) );
CE_f_Emp_P99    = sum( (SS_stats.mu_dist_P99(:)>0) .*SS_stats.mu_dist_Emp(:)  .*CEs_f(:) )/sum( (SS_stats.mu_dist_P99(:)>0) .*SS_stats.mu_dist_Emp(:) );
CE_f_Emp_P999   = sum( (SS_stats.mu_dist_P999(:)>0).*SS_stats.mu_dist_Emp(:)  .*CEs_f(:) )/sum( (SS_stats.mu_dist_P999(:)>0).*SS_stats.mu_dist_Emp(:) );

CE_f_Unemp_Q1   = sum( (SS_stats.mu_dist_Q1(:)>0)  .*SS_stats.mu_dist_Unemp(:)  .*CEs_f(:) )/sum( (SS_stats.mu_dist_Q1(:)>0)  .*SS_stats.mu_dist_Unemp(:) );
CE_f_Unemp_Q2   = sum( (SS_stats.mu_dist_Q2(:)>0)  .*SS_stats.mu_dist_Unemp(:)  .*CEs_f(:) )/sum( (SS_stats.mu_dist_Q2(:)>0)  .*SS_stats.mu_dist_Unemp(:) );
CE_f_Unemp_Q3   = sum( (SS_stats.mu_dist_Q3(:)>0)  .*SS_stats.mu_dist_Unemp(:)  .*CEs_f(:) )/sum( (SS_stats.mu_dist_Q3(:)>0)  .*SS_stats.mu_dist_Unemp(:) );
CE_f_Unemp_Q4   = sum( (SS_stats.mu_dist_Q4(:)>0)  .*SS_stats.mu_dist_Unemp(:)  .*CEs_f(:) )/sum( (SS_stats.mu_dist_Q4(:)>0)  .*SS_stats.mu_dist_Unemp(:) );
CE_f_Unemp_Q5   = sum( (SS_stats.mu_dist_Q5(:)>0)  .*SS_stats.mu_dist_Unemp(:)  .*CEs_f(:) )/sum( (SS_stats.mu_dist_Q5(:)>0)  .*SS_stats.mu_dist_Unemp(:) );
CE_f_Unemp_P90  = sum( (SS_stats.mu_dist_P90(:)>0) .*SS_stats.mu_dist_Unemp(:)  .*CEs_f(:) )/sum( (SS_stats.mu_dist_P90(:)>0) .*SS_stats.mu_dist_Unemp(:) );
CE_f_Unemp_P99  = sum( (SS_stats.mu_dist_P99(:)>0) .*SS_stats.mu_dist_Unemp(:)  .*CEs_f(:) )/sum( (SS_stats.mu_dist_P99(:)>0) .*SS_stats.mu_dist_Unemp(:) );
CE_f_Unemp_P999 = sum( (SS_stats.mu_dist_P999(:)>0).*SS_stats.mu_dist_Unemp(:)  .*CEs_f(:) )/sum( (SS_stats.mu_dist_P999(:)>0).*SS_stats.mu_dist_Unemp(:) );



CE_GOV       = sum( SS_stats.mu_dist(:)      .*CEs_GOV(:) )/sum( SS_stats.mu_dist(:) );
CE_GOV_B01   = sum( SS_stats.mu_dist_B01(:)  .*CEs_GOV(:) )/sum( SS_stats.mu_dist_B01(:) );
CE_GOV_B1    = sum( SS_stats.mu_dist_B1(:)   .*CEs_GOV(:) )/sum( SS_stats.mu_dist_B1(:) );
CE_GOV_B10   = sum( SS_stats.mu_dist_B10(:)  .*CEs_GOV(:) )/sum( SS_stats.mu_dist_B10(:) );
CE_GOV_Q1    = sum( SS_stats.mu_dist_Q1(:)   .*CEs_GOV(:) )/sum( SS_stats.mu_dist_Q1(:) );
CE_GOV_Q2    = sum( SS_stats.mu_dist_Q2(:)   .*CEs_GOV(:) )/sum( SS_stats.mu_dist_Q2(:) );
CE_GOV_Q3    = sum( SS_stats.mu_dist_Q3(:)   .*CEs_GOV(:) )/sum( SS_stats.mu_dist_Q3(:) );
CE_GOV_Q4    = sum( SS_stats.mu_dist_Q4(:)   .*CEs_GOV(:) )/sum( SS_stats.mu_dist_Q4(:) );
CE_GOV_Q5    = sum( SS_stats.mu_dist_Q5(:)   .*CEs_GOV(:) )/sum( SS_stats.mu_dist_Q5(:) );
CE_GOV_P90   = sum( SS_stats.mu_dist_P90(:)  .*CEs_GOV(:) )/sum( SS_stats.mu_dist_P90(:) );
CE_GOV_P99   = sum( SS_stats.mu_dist_P99(:)  .*CEs_GOV(:) )/sum( SS_stats.mu_dist_P99(:) );
CE_GOV_P999  = sum( SS_stats.mu_dist_P999(:) .*CEs_GOV(:) )/sum( SS_stats.mu_dist_P999(:) );
CE_GOV_Ent   = sum( SS_stats.mu_dist_Ent(:)  .*CEs_GOV(:) )/sum( SS_stats.mu_dist_Ent(:) );
CE_GOV_Emp   = sum( SS_stats.mu_dist_Emp(:)  .*CEs_GOV(:) )/sum( SS_stats.mu_dist_Emp(:) );
CE_GOV_Unemp = sum( SS_stats.mu_dist_Unemp(:).*CEs_GOV(:) )/sum( SS_stats.mu_dist_Unemp(:) );

CE_I_GOV_B01   = sum( SS_stats.mu_dist_I_B01(:)  .*CEs_GOV(:) )/sum( SS_stats.mu_dist_I_B01(:) );
CE_I_GOV_B1    = sum( SS_stats.mu_dist_I_B1(:)   .*CEs_GOV(:) )/sum( SS_stats.mu_dist_I_B1(:) );
CE_I_GOV_B10   = sum( SS_stats.mu_dist_I_B10(:)  .*CEs_GOV(:) )/sum( SS_stats.mu_dist_I_B10(:) );
CE_I_GOV_Q1    = sum( SS_stats.mu_dist_I_Q1(:)   .*CEs_GOV(:) )/sum( SS_stats.mu_dist_I_Q1(:) );
CE_I_GOV_Q2    = sum( SS_stats.mu_dist_I_Q2(:)   .*CEs_GOV(:) )/sum( SS_stats.mu_dist_I_Q2(:) );
CE_I_GOV_Q3    = sum( SS_stats.mu_dist_I_Q3(:)   .*CEs_GOV(:) )/sum( SS_stats.mu_dist_I_Q3(:) );
CE_I_GOV_Q4    = sum( SS_stats.mu_dist_I_Q4(:)   .*CEs_GOV(:) )/sum( SS_stats.mu_dist_I_Q4(:) );
CE_I_GOV_Q5    = sum( SS_stats.mu_dist_I_Q5(:)   .*CEs_GOV(:) )/sum( SS_stats.mu_dist_I_Q5(:) );
CE_I_GOV_P90   = sum( SS_stats.mu_dist_I_P90(:)  .*CEs_GOV(:) )/sum( SS_stats.mu_dist_I_P90(:) );
CE_I_GOV_P99   = sum( SS_stats.mu_dist_I_P99(:)  .*CEs_GOV(:) )/sum( SS_stats.mu_dist_I_P99(:) );
CE_I_GOV_P999  = sum( SS_stats.mu_dist_I_P999(:) .*CEs_GOV(:) )/sum( SS_stats.mu_dist_I_P999(:) );


CE_GOV_Ent_Q1     = sum( (SS_stats.mu_dist_Q1(:)>0)  .*SS_stats.mu_dist_Ent(:)  .*CEs_GOV(:) )/sum( (SS_stats.mu_dist_Q1(:)>0)  .*SS_stats.mu_dist_Ent(:) );
CE_GOV_Ent_Q2     = sum( (SS_stats.mu_dist_Q2(:)>0)  .*SS_stats.mu_dist_Ent(:)  .*CEs_GOV(:) )/sum( (SS_stats.mu_dist_Q2(:)>0)  .*SS_stats.mu_dist_Ent(:) );
CE_GOV_Ent_Q3     = sum( (SS_stats.mu_dist_Q3(:)>0)  .*SS_stats.mu_dist_Ent(:)  .*CEs_GOV(:) )/sum( (SS_stats.mu_dist_Q3(:)>0)  .*SS_stats.mu_dist_Ent(:) );
CE_GOV_Ent_Q4     = sum( (SS_stats.mu_dist_Q4(:)>0)  .*SS_stats.mu_dist_Ent(:)  .*CEs_GOV(:) )/sum( (SS_stats.mu_dist_Q4(:)>0)  .*SS_stats.mu_dist_Ent(:) );
CE_GOV_Ent_Q5     = sum( (SS_stats.mu_dist_Q5(:)>0)  .*SS_stats.mu_dist_Ent(:)  .*CEs_GOV(:) )/sum( (SS_stats.mu_dist_Q5(:)>0)  .*SS_stats.mu_dist_Ent(:) );
CE_GOV_Ent_P90    = sum( (SS_stats.mu_dist_P90(:)>0) .*SS_stats.mu_dist_Ent(:)  .*CEs_GOV(:) )/sum( (SS_stats.mu_dist_P90(:)>0) .*SS_stats.mu_dist_Ent(:) );
CE_GOV_Ent_P99    = sum( (SS_stats.mu_dist_P99(:)>0) .*SS_stats.mu_dist_Ent(:)  .*CEs_GOV(:) )/sum( (SS_stats.mu_dist_P99(:)>0) .*SS_stats.mu_dist_Ent(:) );
CE_GOV_Ent_P999   = sum( (SS_stats.mu_dist_P999(:)>0).*SS_stats.mu_dist_Ent(:)  .*CEs_GOV(:) )/sum( (SS_stats.mu_dist_P999(:)>0).*SS_stats.mu_dist_Ent(:) );

CE_GOV_Emp_Q1     = sum( (SS_stats.mu_dist_Q1(:)>0)  .*SS_stats.mu_dist_Emp(:)  .*CEs_GOV(:) )/sum( (SS_stats.mu_dist_Q1(:)>0)  .*SS_stats.mu_dist_Emp(:) );
CE_GOV_Emp_Q2     = sum( (SS_stats.mu_dist_Q2(:)>0)  .*SS_stats.mu_dist_Emp(:)  .*CEs_GOV(:) )/sum( (SS_stats.mu_dist_Q2(:)>0)  .*SS_stats.mu_dist_Emp(:) );
CE_GOV_Emp_Q3     = sum( (SS_stats.mu_dist_Q3(:)>0)  .*SS_stats.mu_dist_Emp(:)  .*CEs_GOV(:) )/sum( (SS_stats.mu_dist_Q3(:)>0)  .*SS_stats.mu_dist_Emp(:) );
CE_GOV_Emp_Q4     = sum( (SS_stats.mu_dist_Q4(:)>0)  .*SS_stats.mu_dist_Emp(:)  .*CEs_GOV(:) )/sum( (SS_stats.mu_dist_Q4(:)>0)  .*SS_stats.mu_dist_Emp(:) );
CE_GOV_Emp_Q5     = sum( (SS_stats.mu_dist_Q5(:)>0)  .*SS_stats.mu_dist_Emp(:)  .*CEs_GOV(:) )/sum( (SS_stats.mu_dist_Q5(:)>0)  .*SS_stats.mu_dist_Emp(:) );
CE_GOV_Emp_P90    = sum( (SS_stats.mu_dist_P90(:)>0) .*SS_stats.mu_dist_Emp(:)  .*CEs_GOV(:) )/sum( (SS_stats.mu_dist_P90(:)>0) .*SS_stats.mu_dist_Emp(:) );
CE_GOV_Emp_P99    = sum( (SS_stats.mu_dist_P99(:)>0) .*SS_stats.mu_dist_Emp(:)  .*CEs_GOV(:) )/sum( (SS_stats.mu_dist_P99(:)>0) .*SS_stats.mu_dist_Emp(:) );
CE_GOV_Emp_P999   = sum( (SS_stats.mu_dist_P999(:)>0).*SS_stats.mu_dist_Emp(:)  .*CEs_GOV(:) )/sum( (SS_stats.mu_dist_P999(:)>0).*SS_stats.mu_dist_Emp(:) );

CE_GOV_Unemp_Q1   = sum( (SS_stats.mu_dist_Q1(:)>0)  .*SS_stats.mu_dist_Unemp(:)  .*CEs_GOV(:) )/sum( (SS_stats.mu_dist_Q1(:)>0)  .*SS_stats.mu_dist_Unemp(:) );
CE_GOV_Unemp_Q2   = sum( (SS_stats.mu_dist_Q2(:)>0)  .*SS_stats.mu_dist_Unemp(:)  .*CEs_GOV(:) )/sum( (SS_stats.mu_dist_Q2(:)>0)  .*SS_stats.mu_dist_Unemp(:) );
CE_GOV_Unemp_Q3   = sum( (SS_stats.mu_dist_Q3(:)>0)  .*SS_stats.mu_dist_Unemp(:)  .*CEs_GOV(:) )/sum( (SS_stats.mu_dist_Q3(:)>0)  .*SS_stats.mu_dist_Unemp(:) );
CE_GOV_Unemp_Q4   = sum( (SS_stats.mu_dist_Q4(:)>0)  .*SS_stats.mu_dist_Unemp(:)  .*CEs_GOV(:) )/sum( (SS_stats.mu_dist_Q4(:)>0)  .*SS_stats.mu_dist_Unemp(:) );
CE_GOV_Unemp_Q5   = sum( (SS_stats.mu_dist_Q5(:)>0)  .*SS_stats.mu_dist_Unemp(:)  .*CEs_GOV(:) )/sum( (SS_stats.mu_dist_Q5(:)>0)  .*SS_stats.mu_dist_Unemp(:) );
CE_GOV_Unemp_P90  = sum( (SS_stats.mu_dist_P90(:)>0) .*SS_stats.mu_dist_Unemp(:)  .*CEs_GOV(:) )/sum( (SS_stats.mu_dist_P90(:)>0) .*SS_stats.mu_dist_Unemp(:) );
CE_GOV_Unemp_P99  = sum( (SS_stats.mu_dist_P99(:)>0) .*SS_stats.mu_dist_Unemp(:)  .*CEs_GOV(:) )/sum( (SS_stats.mu_dist_P99(:)>0) .*SS_stats.mu_dist_Unemp(:) );
CE_GOV_Unemp_P999 = sum( (SS_stats.mu_dist_P999(:)>0).*SS_stats.mu_dist_Unemp(:)  .*CEs_GOV(:) )/sum( (SS_stats.mu_dist_P999(:)>0).*SS_stats.mu_dist_Unemp(:) );

CE_Q       = sum( SS_stats.mu_dist(:)      .*CEs_Q(:) )/sum( SS_stats.mu_dist(:) );
CE_Q_B01   = sum( SS_stats.mu_dist_B01(:)  .*CEs_Q(:) )/sum( SS_stats.mu_dist_B01(:) );
CE_Q_B1    = sum( SS_stats.mu_dist_B1(:)   .*CEs_Q(:) )/sum( SS_stats.mu_dist_B1(:) );
CE_Q_B10   = sum( SS_stats.mu_dist_B10(:)  .*CEs_Q(:) )/sum( SS_stats.mu_dist_B10(:) );
CE_Q_Q1    = sum( SS_stats.mu_dist_Q1(:)   .*CEs_Q(:) )/sum( SS_stats.mu_dist_Q1(:) );
CE_Q_Q2    = sum( SS_stats.mu_dist_Q2(:)   .*CEs_Q(:) )/sum( SS_stats.mu_dist_Q2(:) );
CE_Q_Q3    = sum( SS_stats.mu_dist_Q3(:)   .*CEs_Q(:) )/sum( SS_stats.mu_dist_Q3(:) );
CE_Q_Q4    = sum( SS_stats.mu_dist_Q4(:)   .*CEs_Q(:) )/sum( SS_stats.mu_dist_Q4(:) );
CE_Q_Q5    = sum( SS_stats.mu_dist_Q5(:)   .*CEs_Q(:) )/sum( SS_stats.mu_dist_Q5(:) );
CE_Q_P90   = sum( SS_stats.mu_dist_P90(:)  .*CEs_Q(:) )/sum( SS_stats.mu_dist_P90(:) );
CE_Q_P99   = sum( SS_stats.mu_dist_P99(:)  .*CEs_Q(:) )/sum( SS_stats.mu_dist_P99(:) );
CE_Q_P999  = sum( SS_stats.mu_dist_P999(:) .*CEs_Q(:) )/sum( SS_stats.mu_dist_P999(:) );
CE_Q_Ent   = sum( SS_stats.mu_dist_Ent(:)  .*CEs_Q(:) )/sum( SS_stats.mu_dist_Ent(:) );
CE_Q_Emp   = sum( SS_stats.mu_dist_Emp(:)  .*CEs_Q(:) )/sum( SS_stats.mu_dist_Emp(:) );
CE_Q_Unemp = sum( SS_stats.mu_dist_Unemp(:).*CEs_Q(:) )/sum( SS_stats.mu_dist_Unemp(:) );

CE_I_Q_B01   = sum( SS_stats.mu_dist_I_B01(:)  .*CEs_Q(:) )/sum( SS_stats.mu_dist_I_B01(:) );
CE_I_Q_B1    = sum( SS_stats.mu_dist_I_B1(:)   .*CEs_Q(:) )/sum( SS_stats.mu_dist_I_B1(:) );
CE_I_Q_B10   = sum( SS_stats.mu_dist_I_B10(:)  .*CEs_Q(:) )/sum( SS_stats.mu_dist_I_B10(:) );
CE_I_Q_Q1    = sum( SS_stats.mu_dist_I_Q1(:)   .*CEs_Q(:) )/sum( SS_stats.mu_dist_I_Q1(:) );
CE_I_Q_Q2    = sum( SS_stats.mu_dist_I_Q2(:)   .*CEs_Q(:) )/sum( SS_stats.mu_dist_I_Q2(:) );
CE_I_Q_Q3    = sum( SS_stats.mu_dist_I_Q3(:)   .*CEs_Q(:) )/sum( SS_stats.mu_dist_I_Q3(:) );
CE_I_Q_Q4    = sum( SS_stats.mu_dist_I_Q4(:)   .*CEs_Q(:) )/sum( SS_stats.mu_dist_I_Q4(:) );
CE_I_Q_Q5    = sum( SS_stats.mu_dist_I_Q5(:)   .*CEs_Q(:) )/sum( SS_stats.mu_dist_I_Q5(:) );
CE_I_Q_P90   = sum( SS_stats.mu_dist_I_P90(:)  .*CEs_Q(:) )/sum( SS_stats.mu_dist_I_P90(:) );
CE_I_Q_P99   = sum( SS_stats.mu_dist_I_P99(:)  .*CEs_Q(:) )/sum( SS_stats.mu_dist_I_P99(:) );
CE_I_Q_P999  = sum( SS_stats.mu_dist_I_P999(:) .*CEs_Q(:) )/sum( SS_stats.mu_dist_I_P999(:) );


CE_Q_Ent_Q1     = sum( (SS_stats.mu_dist_Q1(:)>0)  .*SS_stats.mu_dist_Ent(:)  .*CEs_Q(:) )/sum( (SS_stats.mu_dist_Q1(:)>0)  .*SS_stats.mu_dist_Ent(:) );
CE_Q_Ent_Q2     = sum( (SS_stats.mu_dist_Q2(:)>0)  .*SS_stats.mu_dist_Ent(:)  .*CEs_Q(:) )/sum( (SS_stats.mu_dist_Q2(:)>0)  .*SS_stats.mu_dist_Ent(:) );
CE_Q_Ent_Q3     = sum( (SS_stats.mu_dist_Q3(:)>0)  .*SS_stats.mu_dist_Ent(:)  .*CEs_Q(:) )/sum( (SS_stats.mu_dist_Q3(:)>0)  .*SS_stats.mu_dist_Ent(:) );
CE_Q_Ent_Q4     = sum( (SS_stats.mu_dist_Q4(:)>0)  .*SS_stats.mu_dist_Ent(:)  .*CEs_Q(:) )/sum( (SS_stats.mu_dist_Q4(:)>0)  .*SS_stats.mu_dist_Ent(:) );
CE_Q_Ent_Q5     = sum( (SS_stats.mu_dist_Q5(:)>0)  .*SS_stats.mu_dist_Ent(:)  .*CEs_Q(:) )/sum( (SS_stats.mu_dist_Q5(:)>0)  .*SS_stats.mu_dist_Ent(:) );
CE_Q_Ent_P90    = sum( (SS_stats.mu_dist_P90(:)>0) .*SS_stats.mu_dist_Ent(:)  .*CEs_Q(:) )/sum( (SS_stats.mu_dist_P90(:)>0) .*SS_stats.mu_dist_Ent(:) );
CE_Q_Ent_P99    = sum( (SS_stats.mu_dist_P99(:)>0) .*SS_stats.mu_dist_Ent(:)  .*CEs_Q(:) )/sum( (SS_stats.mu_dist_P99(:)>0) .*SS_stats.mu_dist_Ent(:) );
CE_Q_Ent_P999   = sum( (SS_stats.mu_dist_P999(:)>0).*SS_stats.mu_dist_Ent(:)  .*CEs_Q(:) )/sum( (SS_stats.mu_dist_P999(:)>0).*SS_stats.mu_dist_Ent(:) );

CE_Q_Emp_Q1     = sum( (SS_stats.mu_dist_Q1(:)>0)  .*SS_stats.mu_dist_Emp(:)  .*CEs_Q(:) )/sum( (SS_stats.mu_dist_Q1(:)>0)  .*SS_stats.mu_dist_Emp(:) );
CE_Q_Emp_Q2     = sum( (SS_stats.mu_dist_Q2(:)>0)  .*SS_stats.mu_dist_Emp(:)  .*CEs_Q(:) )/sum( (SS_stats.mu_dist_Q2(:)>0)  .*SS_stats.mu_dist_Emp(:) );
CE_Q_Emp_Q3     = sum( (SS_stats.mu_dist_Q3(:)>0)  .*SS_stats.mu_dist_Emp(:)  .*CEs_Q(:) )/sum( (SS_stats.mu_dist_Q3(:)>0)  .*SS_stats.mu_dist_Emp(:) );
CE_Q_Emp_Q4     = sum( (SS_stats.mu_dist_Q4(:)>0)  .*SS_stats.mu_dist_Emp(:)  .*CEs_Q(:) )/sum( (SS_stats.mu_dist_Q4(:)>0)  .*SS_stats.mu_dist_Emp(:) );
CE_Q_Emp_Q5     = sum( (SS_stats.mu_dist_Q5(:)>0)  .*SS_stats.mu_dist_Emp(:)  .*CEs_Q(:) )/sum( (SS_stats.mu_dist_Q5(:)>0)  .*SS_stats.mu_dist_Emp(:) );
CE_Q_Emp_P90    = sum( (SS_stats.mu_dist_P90(:)>0) .*SS_stats.mu_dist_Emp(:)  .*CEs_Q(:) )/sum( (SS_stats.mu_dist_P90(:)>0) .*SS_stats.mu_dist_Emp(:) );
CE_Q_Emp_P99    = sum( (SS_stats.mu_dist_P99(:)>0) .*SS_stats.mu_dist_Emp(:)  .*CEs_Q(:) )/sum( (SS_stats.mu_dist_P99(:)>0) .*SS_stats.mu_dist_Emp(:) );
CE_Q_Emp_P999   = sum( (SS_stats.mu_dist_P999(:)>0).*SS_stats.mu_dist_Emp(:)  .*CEs_Q(:) )/sum( (SS_stats.mu_dist_P999(:)>0).*SS_stats.mu_dist_Emp(:) );

CE_Q_Unemp_Q1   = sum( (SS_stats.mu_dist_Q1(:)>0)  .*SS_stats.mu_dist_Unemp(:)  .*CEs_Q(:) )/sum( (SS_stats.mu_dist_Q1(:)>0)  .*SS_stats.mu_dist_Unemp(:) );
CE_Q_Unemp_Q2   = sum( (SS_stats.mu_dist_Q2(:)>0)  .*SS_stats.mu_dist_Unemp(:)  .*CEs_Q(:) )/sum( (SS_stats.mu_dist_Q2(:)>0)  .*SS_stats.mu_dist_Unemp(:) );
CE_Q_Unemp_Q3   = sum( (SS_stats.mu_dist_Q3(:)>0)  .*SS_stats.mu_dist_Unemp(:)  .*CEs_Q(:) )/sum( (SS_stats.mu_dist_Q3(:)>0)  .*SS_stats.mu_dist_Unemp(:) );
CE_Q_Unemp_Q4   = sum( (SS_stats.mu_dist_Q4(:)>0)  .*SS_stats.mu_dist_Unemp(:)  .*CEs_Q(:) )/sum( (SS_stats.mu_dist_Q4(:)>0)  .*SS_stats.mu_dist_Unemp(:) );
CE_Q_Unemp_Q5   = sum( (SS_stats.mu_dist_Q5(:)>0)  .*SS_stats.mu_dist_Unemp(:)  .*CEs_Q(:) )/sum( (SS_stats.mu_dist_Q5(:)>0)  .*SS_stats.mu_dist_Unemp(:) );
CE_Q_Unemp_P90  = sum( (SS_stats.mu_dist_P90(:)>0) .*SS_stats.mu_dist_Unemp(:)  .*CEs_Q(:) )/sum( (SS_stats.mu_dist_P90(:)>0) .*SS_stats.mu_dist_Unemp(:) );
CE_Q_Unemp_P99  = sum( (SS_stats.mu_dist_P99(:)>0) .*SS_stats.mu_dist_Unemp(:)  .*CEs_Q(:) )/sum( (SS_stats.mu_dist_P99(:)>0) .*SS_stats.mu_dist_Unemp(:) );
CE_Q_Unemp_P999 = sum( (SS_stats.mu_dist_P999(:)>0).*SS_stats.mu_dist_Unemp(:)  .*CEs_Q(:) )/sum( (SS_stats.mu_dist_P999(:)>0).*SS_stats.mu_dist_Unemp(:) );

CE_PI       = sum( SS_stats.mu_dist(:)      .*CEs_PI(:) )/sum( SS_stats.mu_dist(:) );
CE_PI_B01   = sum( SS_stats.mu_dist_B01(:)  .*CEs_PI(:) )/sum( SS_stats.mu_dist_B01(:) );
CE_PI_B1    = sum( SS_stats.mu_dist_B1(:)   .*CEs_PI(:) )/sum( SS_stats.mu_dist_B1(:) );
CE_PI_B10   = sum( SS_stats.mu_dist_B10(:)  .*CEs_PI(:) )/sum( SS_stats.mu_dist_B10(:) );
CE_PI_Q1    = sum( SS_stats.mu_dist_Q1(:)   .*CEs_PI(:) )/sum( SS_stats.mu_dist_Q1(:) );
CE_PI_Q2    = sum( SS_stats.mu_dist_Q2(:)   .*CEs_PI(:) )/sum( SS_stats.mu_dist_Q2(:) );
CE_PI_Q3    = sum( SS_stats.mu_dist_Q3(:)   .*CEs_PI(:) )/sum( SS_stats.mu_dist_Q3(:) );
CE_PI_Q4    = sum( SS_stats.mu_dist_Q4(:)   .*CEs_PI(:) )/sum( SS_stats.mu_dist_Q4(:) );
CE_PI_Q5    = sum( SS_stats.mu_dist_Q5(:)   .*CEs_PI(:) )/sum( SS_stats.mu_dist_Q5(:) );
CE_PI_P90   = sum( SS_stats.mu_dist_P90(:)  .*CEs_PI(:) )/sum( SS_stats.mu_dist_P90(:) );
CE_PI_P99   = sum( SS_stats.mu_dist_P99(:)  .*CEs_PI(:) )/sum( SS_stats.mu_dist_P99(:) );
CE_PI_P999  = sum( SS_stats.mu_dist_P999(:) .*CEs_PI(:) )/sum( SS_stats.mu_dist_P999(:) );
CE_PI_Ent   = sum( SS_stats.mu_dist_Ent(:)  .*CEs_PI(:) )/sum( SS_stats.mu_dist_Ent(:) );
CE_PI_Emp   = sum( SS_stats.mu_dist_Emp(:)  .*CEs_PI(:) )/sum( SS_stats.mu_dist_Emp(:) );
CE_PI_Unemp = sum( SS_stats.mu_dist_Unemp(:).*CEs_PI(:) )/sum( SS_stats.mu_dist_Unemp(:) );


CE_PI_W       = sum( SS_stats.mu_dist(:)      .*CEs_PI_W(:) )/sum( SS_stats.mu_dist(:) );
CE_PI_W_B01   = sum( SS_stats.mu_dist_B01(:)  .*CEs_PI_W(:) )/sum( SS_stats.mu_dist_B01(:) );
CE_PI_W_B1    = sum( SS_stats.mu_dist_B1(:)   .*CEs_PI_W(:) )/sum( SS_stats.mu_dist_B1(:) );
CE_PI_W_B10   = sum( SS_stats.mu_dist_B10(:)  .*CEs_PI_W(:) )/sum( SS_stats.mu_dist_B10(:) );
CE_PI_W_Q1    = sum( SS_stats.mu_dist_Q1(:)   .*CEs_PI_W(:) )/sum( SS_stats.mu_dist_Q1(:) );
CE_PI_W_Q2    = sum( SS_stats.mu_dist_Q2(:)   .*CEs_PI_W(:) )/sum( SS_stats.mu_dist_Q2(:) );
CE_PI_W_Q3    = sum( SS_stats.mu_dist_Q3(:)   .*CEs_PI_W(:) )/sum( SS_stats.mu_dist_Q3(:) );
CE_PI_W_Q4    = sum( SS_stats.mu_dist_Q4(:)   .*CEs_PI_W(:) )/sum( SS_stats.mu_dist_Q4(:) );
CE_PI_W_Q5    = sum( SS_stats.mu_dist_Q5(:)   .*CEs_PI_W(:) )/sum( SS_stats.mu_dist_Q5(:) );
CE_PI_W_P90   = sum( SS_stats.mu_dist_P90(:)  .*CEs_PI_W(:) )/sum( SS_stats.mu_dist_P90(:) );
CE_PI_W_P99   = sum( SS_stats.mu_dist_P99(:)  .*CEs_PI_W(:) )/sum( SS_stats.mu_dist_P99(:) );
CE_PI_W_P999  = sum( SS_stats.mu_dist_P999(:) .*CEs_PI_W(:) )/sum( SS_stats.mu_dist_P999(:) );
CE_PI_W_Ent   = sum( SS_stats.mu_dist_Ent(:)  .*CEs_PI_W(:) )/sum( SS_stats.mu_dist_Ent(:) );
CE_PI_W_Emp   = sum( SS_stats.mu_dist_Emp(:)  .*CEs_PI_W(:) )/sum( SS_stats.mu_dist_Emp(:) );
CE_PI_W_Unemp = sum( SS_stats.mu_dist_Unemp(:).*CEs_PI_W(:) )/sum( SS_stats.mu_dist_Unemp(:) );

CE_PI_R       = sum( SS_stats.mu_dist(:)      .*CEs_PI_R(:) )/sum( SS_stats.mu_dist(:) );
CE_PI_R_B01   = sum( SS_stats.mu_dist_B01(:)  .*CEs_PI_R(:) )/sum( SS_stats.mu_dist_B01(:) );
CE_PI_R_B1    = sum( SS_stats.mu_dist_B1(:)   .*CEs_PI_R(:) )/sum( SS_stats.mu_dist_B1(:) );
CE_PI_R_B10   = sum( SS_stats.mu_dist_B10(:)  .*CEs_PI_R(:) )/sum( SS_stats.mu_dist_B10(:) );
CE_PI_R_Q1    = sum( SS_stats.mu_dist_Q1(:)   .*CEs_PI_R(:) )/sum( SS_stats.mu_dist_Q1(:) );
CE_PI_R_Q2    = sum( SS_stats.mu_dist_Q2(:)   .*CEs_PI_R(:) )/sum( SS_stats.mu_dist_Q2(:) );
CE_PI_R_Q3    = sum( SS_stats.mu_dist_Q3(:)   .*CEs_PI_R(:) )/sum( SS_stats.mu_dist_Q3(:) );
CE_PI_R_Q4    = sum( SS_stats.mu_dist_Q4(:)   .*CEs_PI_R(:) )/sum( SS_stats.mu_dist_Q4(:) );
CE_PI_R_Q5    = sum( SS_stats.mu_dist_Q5(:)   .*CEs_PI_R(:) )/sum( SS_stats.mu_dist_Q5(:) );
CE_PI_R_P90   = sum( SS_stats.mu_dist_P90(:)  .*CEs_PI_R(:) )/sum( SS_stats.mu_dist_P90(:) );
CE_PI_R_P99   = sum( SS_stats.mu_dist_P99(:)  .*CEs_PI_R(:) )/sum( SS_stats.mu_dist_P99(:) );
CE_PI_R_P999  = sum( SS_stats.mu_dist_P999(:) .*CEs_PI_R(:) )/sum( SS_stats.mu_dist_P999(:) );
CE_PI_R_Ent   = sum( SS_stats.mu_dist_Ent(:)  .*CEs_PI_R(:) )/sum( SS_stats.mu_dist_Ent(:) );
CE_PI_R_Emp   = sum( SS_stats.mu_dist_Emp(:)  .*CEs_PI_R(:) )/sum( SS_stats.mu_dist_Emp(:) );
CE_PI_R_Unemp = sum( SS_stats.mu_dist_Unemp(:).*CEs_PI_R(:) )/sum( SS_stats.mu_dist_Unemp(:) );

CE_I_PI_B01   = sum( SS_stats.mu_dist_I_B01(:)  .*CEs_PI(:) )/sum( SS_stats.mu_dist_I_B01(:) );
CE_I_PI_B1    = sum( SS_stats.mu_dist_I_B1(:)   .*CEs_PI(:) )/sum( SS_stats.mu_dist_I_B1(:) );
CE_I_PI_B10   = sum( SS_stats.mu_dist_I_B10(:)  .*CEs_PI(:) )/sum( SS_stats.mu_dist_I_B10(:) );
CE_I_PI_Q1    = sum( SS_stats.mu_dist_I_Q1(:)   .*CEs_PI(:) )/sum( SS_stats.mu_dist_I_Q1(:) );
CE_I_PI_Q2    = sum( SS_stats.mu_dist_I_Q2(:)   .*CEs_PI(:) )/sum( SS_stats.mu_dist_I_Q2(:) );
CE_I_PI_Q3    = sum( SS_stats.mu_dist_I_Q3(:)   .*CEs_PI(:) )/sum( SS_stats.mu_dist_I_Q3(:) );
CE_I_PI_Q4    = sum( SS_stats.mu_dist_I_Q4(:)   .*CEs_PI(:) )/sum( SS_stats.mu_dist_I_Q4(:) );
CE_I_PI_Q5    = sum( SS_stats.mu_dist_I_Q5(:)   .*CEs_PI(:) )/sum( SS_stats.mu_dist_I_Q5(:) );
CE_I_PI_P90   = sum( SS_stats.mu_dist_I_P90(:)  .*CEs_PI(:) )/sum( SS_stats.mu_dist_I_P90(:) );
CE_I_PI_P99   = sum( SS_stats.mu_dist_I_P99(:)  .*CEs_PI(:) )/sum( SS_stats.mu_dist_I_P99(:) );
CE_I_PI_P999  = sum( SS_stats.mu_dist_I_P999(:) .*CEs_PI(:) )/sum( SS_stats.mu_dist_I_P999(:) );

CE_I_PI_W_B01   = sum( SS_stats.mu_dist_I_B01(:)  .*CEs_PI_W(:) )/sum( SS_stats.mu_dist_I_B01(:) );
CE_I_PI_W_B1    = sum( SS_stats.mu_dist_I_B1(:)   .*CEs_PI_W(:) )/sum( SS_stats.mu_dist_I_B1(:) );
CE_I_PI_W_B10   = sum( SS_stats.mu_dist_I_B10(:)  .*CEs_PI_W(:) )/sum( SS_stats.mu_dist_I_B10(:) );
CE_I_PI_W_Q1    = sum( SS_stats.mu_dist_I_Q1(:)   .*CEs_PI_W(:) )/sum( SS_stats.mu_dist_I_Q1(:) );
CE_I_PI_W_Q2    = sum( SS_stats.mu_dist_I_Q2(:)   .*CEs_PI_W(:) )/sum( SS_stats.mu_dist_I_Q2(:) );
CE_I_PI_W_Q3    = sum( SS_stats.mu_dist_I_Q3(:)   .*CEs_PI_W(:) )/sum( SS_stats.mu_dist_I_Q3(:) );
CE_I_PI_W_Q4    = sum( SS_stats.mu_dist_I_Q4(:)   .*CEs_PI_W(:) )/sum( SS_stats.mu_dist_I_Q4(:) );
CE_I_PI_W_Q5    = sum( SS_stats.mu_dist_I_Q5(:)   .*CEs_PI_W(:) )/sum( SS_stats.mu_dist_I_Q5(:) );
CE_I_PI_W_P90   = sum( SS_stats.mu_dist_I_P90(:)  .*CEs_PI_W(:) )/sum( SS_stats.mu_dist_I_P90(:) );
CE_I_PI_W_P99   = sum( SS_stats.mu_dist_I_P99(:)  .*CEs_PI_W(:) )/sum( SS_stats.mu_dist_I_P99(:) );
CE_I_PI_W_P999  = sum( SS_stats.mu_dist_I_P999(:) .*CEs_PI_W(:) )/sum( SS_stats.mu_dist_I_P999(:) );

CE_I_PI_R_B01   = sum( SS_stats.mu_dist_I_B01(:)  .*CEs_PI_R(:) )/sum( SS_stats.mu_dist_I_B01(:) );
CE_I_PI_R_B1    = sum( SS_stats.mu_dist_I_B1(:)   .*CEs_PI_R(:) )/sum( SS_stats.mu_dist_I_B1(:) );
CE_I_PI_R_B10   = sum( SS_stats.mu_dist_I_B10(:)  .*CEs_PI_R(:) )/sum( SS_stats.mu_dist_I_B10(:) );
CE_I_PI_R_Q1    = sum( SS_stats.mu_dist_I_Q1(:)   .*CEs_PI_R(:) )/sum( SS_stats.mu_dist_I_Q1(:) );
CE_I_PI_R_Q2    = sum( SS_stats.mu_dist_I_Q2(:)   .*CEs_PI_R(:) )/sum( SS_stats.mu_dist_I_Q2(:) );
CE_I_PI_R_Q3    = sum( SS_stats.mu_dist_I_Q3(:)   .*CEs_PI_R(:) )/sum( SS_stats.mu_dist_I_Q3(:) );
CE_I_PI_R_Q4    = sum( SS_stats.mu_dist_I_Q4(:)   .*CEs_PI_R(:) )/sum( SS_stats.mu_dist_I_Q4(:) );
CE_I_PI_R_Q5    = sum( SS_stats.mu_dist_I_Q5(:)   .*CEs_PI_R(:) )/sum( SS_stats.mu_dist_I_Q5(:) );
CE_I_PI_R_P90   = sum( SS_stats.mu_dist_I_P90(:)  .*CEs_PI_R(:) )/sum( SS_stats.mu_dist_I_P90(:) );
CE_I_PI_R_P99   = sum( SS_stats.mu_dist_I_P99(:)  .*CEs_PI_R(:) )/sum( SS_stats.mu_dist_I_P99(:) );
CE_I_PI_R_P999  = sum( SS_stats.mu_dist_I_P999(:) .*CEs_PI_R(:) )/sum( SS_stats.mu_dist_I_P999(:) );

CE_PI_collect   = [CE_PI_B01,CE_PI_B1,CE_PI_B10,CE_PI_Q1,CE_PI_Q2,CE_PI_Q3,CE_PI_Q4,CE_PI_Q5,CE_PI_P90,CE_PI_P99,CE_PI_P999];
CE_PI_W_collect = [CE_PI_W_B01,CE_PI_W_B1,CE_PI_W_B10,CE_PI_W_Q1,CE_PI_W_Q2,CE_PI_W_Q3,CE_PI_W_Q4,CE_PI_W_Q5,CE_PI_W_P90,CE_PI_W_P99,CE_PI_W_P999];
CE_PI_R_collect = [CE_PI_R_B01,CE_PI_R_B1,CE_PI_R_B10,CE_PI_R_Q1,CE_PI_R_Q2,CE_PI_R_Q3,CE_PI_R_Q4,CE_PI_R_Q5,CE_PI_R_P90,CE_PI_R_P99,CE_PI_R_P999];

CE_I_PI_collect   = [CE_I_PI_B01,CE_I_PI_B1,CE_I_PI_B10,CE_I_PI_Q1,CE_I_PI_Q2,CE_I_PI_Q3,CE_I_PI_Q4,CE_I_PI_Q5,CE_I_PI_P90,CE_I_PI_P99,CE_I_PI_P999];
CE_I_PI_W_collect = [CE_I_PI_W_B01,CE_I_PI_W_B1,CE_I_PI_W_B10,CE_I_PI_W_Q1,CE_I_PI_W_Q2,CE_I_PI_W_Q3,CE_I_PI_W_Q4,CE_I_PI_W_Q5,CE_I_PI_W_P90,CE_I_PI_W_P99,CE_I_PI_W_P999];
CE_I_PI_R_collect = [CE_I_PI_R_B01,CE_I_PI_R_B1,CE_I_PI_R_B10,CE_I_PI_R_Q1,CE_I_PI_R_Q2,CE_I_PI_R_Q3,CE_I_PI_R_Q4,CE_I_PI_R_Q5,CE_I_PI_R_P90,CE_I_PI_R_P99,CE_I_PI_R_P999];

CE_PI_bar = [CE_PI_collect;CE_PI_W_collect;CE_PI_R_collect];

CE_I_PI_bar = [CE_I_PI_collect;CE_I_PI_W_collect;CE_I_PI_R_collect];

CE_collect          = [CE_B01,CE_B1,CE_B10,CE_Q1,CE_Q2,CE_Q3,CE_Q4,CE_Q5,CE_P90,CE_P99,CE_P999];
CE_RR_collect       = [CE_RR_B01,CE_RR_B1,CE_RR_B10,CE_RR_Q1,CE_RR_Q2,CE_RR_Q3,CE_RR_Q4,CE_RR_Q5,CE_RR_P90,CE_RR_P99,CE_RR_P999];
CE_PROFIT_collect   = [CE_PROFIT_B01,CE_PROFIT_B1,CE_PROFIT_B10,CE_PROFIT_Q1,CE_PROFIT_Q2,CE_PROFIT_Q3,CE_PROFIT_Q4,CE_PROFIT_Q5,CE_PROFIT_P90,CE_PROFIT_P99,CE_PROFIT_P999];
CE_Q_collect        = [CE_Q_B01,CE_Q_B1,CE_Q_B10,CE_Q_Q1,CE_Q_Q2,CE_Q_Q3,CE_Q_Q4,CE_Q_Q5,CE_Q_P90,CE_Q_P99,CE_Q_P999];
CE_W_collect        = [CE_W_B01,CE_W_B1,CE_W_B10,CE_W_Q1,CE_W_Q2,CE_W_Q3,CE_W_Q4,CE_W_Q5,CE_W_P90,CE_W_P99,CE_W_P999];
CE_f_collect        = [CE_f_B01,CE_f_B1,CE_f_B10,CE_f_Q1,CE_f_Q2,CE_f_Q3,CE_f_Q4,CE_f_Q5,CE_f_P90,CE_f_P99,CE_f_P999];
CE_GOV_collect      = [CE_GOV_B01,CE_GOV_B1,CE_GOV_B10,CE_GOV_Q1,CE_GOV_Q2,CE_GOV_Q3,CE_GOV_Q4,CE_GOV_Q5,CE_GOV_P90,CE_GOV_P99,CE_GOV_P999];

CE_bar = [CE_collect;CE_RR_collect;CE_PROFIT_collect;CE_Q_collect;CE_W_collect;CE_f_collect;CE_GOV_collect];

CE_I_collect          = [CE_I_B01,CE_I_B1,CE_I_B10,CE_I_Q1,CE_I_Q2,CE_I_Q3,CE_I_Q4,CE_I_Q5,CE_I_P90,CE_I_P99,CE_I_P999];
CE_I_RR_collect       = [CE_I_RR_B01,CE_I_RR_B1,CE_I_RR_B10,CE_I_RR_Q1,CE_I_RR_Q2,CE_I_RR_Q3,CE_I_RR_Q4,CE_I_RR_Q5,CE_I_RR_P90,CE_I_RR_P99,CE_I_RR_P999];
CE_I_PROFIT_collect   = [CE_I_PROFIT_B01,CE_I_PROFIT_B1,CE_I_PROFIT_B10,CE_I_PROFIT_Q1,CE_I_PROFIT_Q2,CE_I_PROFIT_Q3,CE_I_PROFIT_Q4,CE_I_PROFIT_Q5,CE_I_PROFIT_P90,CE_I_PROFIT_P99,CE_I_PROFIT_P999];
CE_I_Q_collect        = [CE_I_Q_B01,CE_I_Q_B1,CE_I_Q_B10,CE_I_Q_Q1,CE_I_Q_Q2,CE_I_Q_Q3,CE_I_Q_Q4,CE_I_Q_Q5,CE_I_Q_P90,CE_I_Q_P99,CE_I_Q_P999];
CE_I_W_collect        = [CE_I_W_B01,CE_I_W_B1,CE_I_W_B10,CE_I_W_Q1,CE_I_W_Q2,CE_I_W_Q3,CE_I_W_Q4,CE_I_W_Q5,CE_I_W_P90,CE_I_W_P99,CE_I_W_P999];
CE_I_f_collect        = [CE_I_f_B01,CE_I_f_B1,CE_I_f_B10,CE_I_f_Q1,CE_I_f_Q2,CE_I_f_Q3,CE_I_f_Q4,CE_I_f_Q5,CE_I_f_P90,CE_I_f_P99,CE_I_f_P999];
CE_I_GOV_collect      = [CE_I_GOV_B01,CE_I_GOV_B1,CE_I_GOV_B10,CE_I_GOV_Q1,CE_I_GOV_Q2,CE_I_GOV_Q3,CE_I_GOV_Q4,CE_I_GOV_Q5,CE_I_GOV_P90,CE_I_GOV_P99,CE_I_GOV_P999];

CE_I_bar = [CE_I_collect;CE_I_RR_collect;CE_I_PROFIT_collect;CE_I_Q_collect;CE_I_W_collect;CE_I_f_collect;CE_I_GOV_collect];

