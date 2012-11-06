/* @HEADER@ */
// ************************************************************************
//
//                              Fitob
//            Copyright (2012) Janos Benk (benkjanos@gmail.com)
//
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Janos Benk (benkjanos@gmail.com),
//
// ************************************************************************
/* @HEADER@ */

/** This files defines the scripting grammar */

#ifndef FITOB_SCRIPT_GRAMMAR_HPP_
#define FITOB_SCRIPT_GRAMMAR_HPP_

using namespace BOOST_SPIRIT_CLASSIC_NS;

namespace fitob{

  /** Here we specify the grammar of the Theta scipt
   * The result will be one XML tree which can later be processed
   * and the according object structure will be created */
  struct ThetaParser : public grammar<ThetaParser>
  {
    static const int constantID = 1;
    static const int factorID = 2;
    static const int termID = 3;
    static const int expressionID = 4;
    static const int variableID = 5;
    static const int endtermID = 6;
    static const int minFctID = 7;
    static const int maxFctID = 8;
    static const int sqrtFctID = 9;
    static const int expFctID = 10;
    static const int logFctID = 11;
    static const int expectedFctID = 12;
    static const int discountFctID = 13;
    static const int condFctID = 14;
    
    static const int BodyID = 30;
    static const int OpID = 31;
    static const int AssigmentOpID = 32;
    static const int ThetaOpID = 33;
    static const int CommentOpID = 34;
    static const int LoopOpID = 35;
    static const int IfElseOpID = 36;
    static const int IfOpID = 37;
    static const int PlOpID = 38;
    
    static const int ScriptID = 50;
    static const int DeclarationID = 51;
    
    static const int SDEDeclarationID = 60;
    static const int SDEModelDefID = 61;

    static const int modesimplID = 65;
    static const int modefullID = 66;
    static const int modeInterestID = 67;

    template <typename ScannerT>
    struct definition
    {
        definition(ThetaParser const& /*self*/)
        {
            //  Start grammar definition
	    
	    // ---------- Define Expressions ------------
            endterm     =   constant
	                    |   minFct
	                    |   maxFct
                        |   sqrtFct
                        |   expFct
                        |   logFct
                        |   expectedFct
                        |   discountFct
                        |   condFct
	                    |   variable
                        |   inner_node_d[ch_p('(') >> expression >> ch_p(')')]
                        |   (root_node_d[ch_p('-')] >> endterm);

            constant     =   leaf_node_d[ 
                               (!ch_p('-') >> discard_node_d[*(blank_p)] 
                               >> lexeme_d[real_p] 
                               >> discard_node_d[*(blank_p)])
                            ];

            variable    =   leaf_node_d[
                              (!ch_p('-') >> discard_node_d[*(blank_p)] 
                              >> lexeme_d[+(alnum_p) >> !ch_p('!')] >>
                              discard_node_d[*(blank_p)])
                            ];

            factor      =   endterm >>
                            *(  (root_node_d[ch_p('^')] >> endterm)
                            );

            term        =   factor >>
                            *(  (root_node_d[ch_p('*')] >> factor)
                              | (root_node_d[ch_p('/')] >> factor)
                            );

            expression  =   (term >>
                               *(  
                                     (root_node_d[lexeme_d[chseq_p("<=")]] >> term)
                                   | (root_node_d[lexeme_d[chseq_p(">=")]] >> term)
                                   | (root_node_d[lexeme_d[chseq_p("==")]] >> term)
                                   | (root_node_d[lexeme_d[chseq_p("~=")]] >> term)
                                   | (root_node_d[lexeme_d[chseq_p("||")]] >> term)
                                   | (root_node_d[lexeme_d[chseq_p("&&")]] >> term)
                                   | (root_node_d[ch_p('>')] >> term)
                                   | (root_node_d[ch_p('<')] >> term)
                                   | (root_node_d[ch_p('+')] >> term)
                                   | (root_node_d[ch_p('-')] >> term)
                                )
                            );

            minFct      =  (discard_node_d[*(blank_p)] >> lexeme_d[  //infix deletse the , from the middle, inner_node_d deletese the brackets
                             ( root_node_d[chseq_p("MIN")] >> 
                               infix_node_d[ inner_node_d[ ch_p('(') >> expression >>
                                           ch_p(',') >> expression >> ch_p(')') ]])
                               ]
                            );
    
            maxFct      =  (discard_node_d[*(blank_p)] >> lexeme_d[  //infix deletse the , from the middle, inner_node_d deletese the brackets
                             ( root_node_d[chseq_p("MAX")] >> 
                               infix_node_d[ inner_node_d[ ch_p('(') >> expression >>
                                           ch_p(',') >> expression >> ch_p(')') ]])
                               ]
                            );

            sqrtFct     =  (discard_node_d[*(blank_p)] >> lexeme_d[  //infix deletse the , from the middle, inner_node_d deletese the brackets
                             ( root_node_d[chseq_p("SQRT")] >> 
                               inner_node_d[ ch_p('(') >> expression >> ch_p(')') ])
                               ]
                            );

            expFct     =  (discard_node_d[*(blank_p)] >> lexeme_d[  //infix deletse the , from the middle, inner_node_d deletese the brackets
                             ( root_node_d[chseq_p("EXP")] >> 
                               inner_node_d[ ch_p('(') >> expression >> ch_p(')') ])
                               ]
                            );

            logFct     =  (discard_node_d[*(blank_p)] >> lexeme_d[  //infix deletse the , from the middle, inner_node_d deletese the brackets
                             ( root_node_d[chseq_p("LOG")] >> 
                               inner_node_d[ ch_p('(') >> expression >> ch_p(')') ])
                               ]
                            );

            expectedFct  =  (discard_node_d[*(blank_p)] >> lexeme_d[  //infix deletse the , from the middle, inner_node_d deletese the brackets
                             ( root_node_d[chseq_p("EXPECT")] >>
                               inner_node_d[ ch_p('(') >> expression >> ch_p(')') ])
                               ]
                            );

            discountFct  =  (discard_node_d[*(blank_p)] >> lexeme_d[  //infix deletse the , from the middle, inner_node_d deletese the brackets
                              ( root_node_d[chseq_p("DISCOUNT")] >>
                                infix_node_d[ inner_node_d[ ch_p('(') >> expression >>
                                            ch_p(',') >> expression >> ch_p(')') ]])
                                ]
                             );

            condFct  =  (discard_node_d[*(blank_p)] >> lexeme_d[  //infix deletse the , from the middle, inner_node_d deletese the brackets
                              ( root_node_d[chseq_p("COND")] >>
                                infix_node_d[ inner_node_d[ ch_p('(') >> expression >>
                                            ch_p(',') >> expression >> ch_p(',') >> expression >> ch_p(')') ]])
                                ]
                             );
	    // ---------- END Define Expressions ------------
	    
        // ---------- Operators definition --------------

	    Body        = *(Op);
	    
	    Op          =  discard_node_d[*(blank_p)] >> 
                           (  ThetaOp
                            | LoopOp 
                            | IfElseOp
                            | IfOp
                            | PlOp
                            | AssigmentOp
                            | CommentOp
                            ) >> 
                            discard_node_d[*(blank_p)] >> 
                            discard_node_d[( +(eol_p) | +(ch_p(';')))]
                            >> discard_node_d[*(blank_p)]
                            //>> discard_node_d[ +(ch_p('%')) >> *(alnum_p) >> +(eol_p) ] 
                            ;

	    CommentOp   =  discard_node_d[*(blank_p) >> +(ch_p('%')) >> +(alnum_p)];
	    
	    AssigmentOp =  ( variable >> root_node_d[ch_p('=')] >>
	    		         discard_node_d[*(blank_p)] >> expression);

	        PlOp   =  ( root_node_d[lexeme_d[chseq_p("Plot")]] >> discard_node_d[+(blank_p)] >> expression );

            ThetaOp     =  (root_node_d[lexeme_d[chseq_p("Theta")]] >> discard_node_d[+(blank_p)] >> expression);
	    
            LoopOp     =  (root_node_d[lexeme_d[chseq_p("loop")]] >> discard_node_d[+(blank_p)] 
                           >> expression >> discard_node_d[+(blank_p)] >> discard_node_d[*(eol_p)] 
                           >> Body 
                           >> discard_node_d[lexeme_d[chseq_p("end")]]);

            IfOp       =  (root_node_d[lexeme_d[chseq_p("if")]] >> discard_node_d[+(blank_p)] 
                           >> expression >> discard_node_d[+(blank_p)] >> discard_node_d[*(eol_p)] 
                           >> Body 
                           >> discard_node_d[lexeme_d[chseq_p("end")]]);

            IfElseOp   =  (root_node_d[lexeme_d[chseq_p("if")]] >> discard_node_d[+(blank_p)] 
                           >> expression >> discard_node_d[+(blank_p)] >> discard_node_d[*(eol_p)] 
                           >> Body 
                           >> discard_node_d[lexeme_d[chseq_p("else")]] >> discard_node_d[+(blank_p)]
                           >> discard_node_d[*(eol_p)] 
                           >> Body
                           >> discard_node_d[lexeme_d[chseq_p("end")]]
                           );

            // ---------- END Operators definition --------------

            // --------------- Script definition ----------------

            Script     =  ( root_node_d[lexeme_d[chseq_p("model")]] 
                            >> discard_node_d[+(blank_p)]
                            >> leaf_node_d[lexeme_d[+(alnum_p)]] >> discard_node_d[*(blank_p)]
                            >> discard_node_d[( +(eol_p) | +(ch_p(';')))]
                            >> +(Declaration)
                            >> *(SDEDeclaration)
                            >> Body
                            >> discard_node_d[lexeme_d[chseq_p("end")]] 
                            >> discard_node_d[(+(eol_p) | +(ch_p(';')))]
                           );

           SDEDeclaration =  (
        		            discard_node_d[*(blank_p)] >> variable >> discard_node_d[*(blank_p)]
        		            >> root_node_d[ch_p('=')] >> discard_node_d[*(blank_p)]
        		            >> modefull
        		 		     | modeInterest
        		 		     | modesimpl
        		            >> discard_node_d[( +(eol_p) | +(ch_p(';')))]
        		           );
           // todo: this is not complete yet, work out how to define a general interest rate model
           //SDEModelDef =    modesimpl
		   //                | modeInterest
		   //                | modefull ;

           modefull    =  (  //infix deletse the , from the middle, inner_node_d deletese the brackets
                             ( lexeme_d[root_node_d[chseq_p("MODEL")]] >>
                               infix_node_d[ inner_node_d[ ch_p('(') >> expression >>
                               ch_p(',') >> expression >> ch_p(',') >> expression >> ch_p(')') ]] >> discard_node_d[ch_p(';')] )
                           );

           modesimpl    =  (  //infix deletse the , from the middle, inner_node_d deletese the brackets
                             ( lexeme_d[root_node_d[chseq_p("MODEL")]] >>
                             infix_node_d[ inner_node_d[ ch_p('(') >> expression >>
                             ch_p(',') >> expression >> ch_p(')')] >> discard_node_d[ch_p(';')] ])
                            );

           modeInterest =  ( //infix deletse the , from the middle, inner_node_d deletese the brackets
                           ( lexeme_d[root_node_d[chseq_p("INTEREST")]] >>
                            infix_node_d[ inner_node_d[ ch_p('(') >> expression >>
                            ch_p(',') >> expression >> ch_p(',') >> expression >> ch_p(')') ]] >> discard_node_d[ch_p(';')] )
                           );
                         
           Declaration =  ( 
                            discard_node_d[*(blank_p)]
                            >> root_node_d[lexeme_d[ chseq_p("import") | chseq_p("export")]]
                            >> discard_node_d[+(blank_p)]
                            >> leaf_node_d[lexeme_d[+(alnum_p)]]
                            >> discard_node_d[( +(eol_p) | +(ch_p(';')))]
	                   );
                                            
            // -----------END Script definition ----------------

            // --------- End grammar definition ------------

            // turn on the debugging info.
            BOOST_SPIRIT_DEBUG_RULE(constant);
            BOOST_SPIRIT_DEBUG_RULE(factor);
            BOOST_SPIRIT_DEBUG_RULE(term);
            BOOST_SPIRIT_DEBUG_RULE(expression);
            BOOST_SPIRIT_DEBUG_RULE(variable);
            BOOST_SPIRIT_DEBUG_RULE(minFct);
            BOOST_SPIRIT_DEBUG_RULE(maxFct);
            BOOST_SPIRIT_DEBUG_RULE(sqrtFct);
            BOOST_SPIRIT_DEBUG_RULE(expFct);
            BOOST_SPIRIT_DEBUG_RULE(logFct);
            BOOST_SPIRIT_DEBUG_RULE(expectedFct);
            BOOST_SPIRIT_DEBUG_RULE(discountFct);
            BOOST_SPIRIT_DEBUG_RULE(condFct);
	    
            BOOST_SPIRIT_DEBUG_RULE(Body);
            BOOST_SPIRIT_DEBUG_RULE(Op);
            BOOST_SPIRIT_DEBUG_RULE(AssigmentOp);
            BOOST_SPIRIT_DEBUG_RULE(ThetaOp);
            BOOST_SPIRIT_DEBUG_RULE(LoopOp);
            BOOST_SPIRIT_DEBUG_RULE(IfOp);
            BOOST_SPIRIT_DEBUG_RULE(IfElseOp);
            BOOST_SPIRIT_DEBUG_RULE(PlOp);
	    
            BOOST_SPIRIT_DEBUG_RULE(Script);
            BOOST_SPIRIT_DEBUG_RULE(Declaration);

            BOOST_SPIRIT_DEBUG_RULE(SDEDeclaration);
            BOOST_SPIRIT_DEBUG_RULE(SDEModelDef);
            BOOST_SPIRIT_DEBUG_RULE(modefull);
            BOOST_SPIRIT_DEBUG_RULE(modesimpl);
            BOOST_SPIRIT_DEBUG_RULE(modeInterest);
        }

        rule<ScannerT, parser_context<>, parser_tag<expressionID> >   expression;
        rule<ScannerT, parser_context<>, parser_tag<termID> >         term;
        rule<ScannerT, parser_context<>, parser_tag<factorID> >       factor;
        rule<ScannerT, parser_context<>, parser_tag<constantID> >     constant;
        rule<ScannerT, parser_context<>, parser_tag<variableID> >     variable;
        rule<ScannerT, parser_context<>, parser_tag<endtermID> >      endterm;
        rule<ScannerT, parser_context<>, parser_tag<minFctID> >       minFct;
        rule<ScannerT, parser_context<>, parser_tag<maxFctID> >       maxFct;
        rule<ScannerT, parser_context<>, parser_tag<sqrtFctID> >      sqrtFct;
        rule<ScannerT, parser_context<>, parser_tag<expFctID> >       expFct;
        rule<ScannerT, parser_context<>, parser_tag<logFctID> >       logFct;
        rule<ScannerT, parser_context<>, parser_tag<expectedFctID> >  expectedFct;
        rule<ScannerT, parser_context<>, parser_tag<discountFctID> >  discountFct;
        rule<ScannerT, parser_context<>, parser_tag<condFctID> >      condFct;
        /* */
        rule<ScannerT, parser_context<>, parser_tag<BodyID> >          Body;
        rule<ScannerT, parser_context<>, parser_tag<OpID> >            Op;
        rule<ScannerT, parser_context<>, parser_tag<AssigmentOpID> >   AssigmentOp;
        rule<ScannerT, parser_context<>, parser_tag<ThetaOpID> >       ThetaOp;
        rule<ScannerT, parser_context<>, parser_tag<CommentOpID> >     CommentOp;
        rule<ScannerT, parser_context<>, parser_tag<LoopOpID> >        LoopOp;
        rule<ScannerT, parser_context<>, parser_tag<IfOpID> >          IfOp;
        rule<ScannerT, parser_context<>, parser_tag<IfElseOpID> >      IfElseOp;
        rule<ScannerT, parser_context<>, parser_tag<PlOpID> >          PlOp;

        rule<ScannerT, parser_context<>, parser_tag<ScriptID> >        Script;
        rule<ScannerT, parser_context<>, parser_tag<DeclarationID> >   Declaration;

        rule<ScannerT, parser_context<>, parser_tag<SDEDeclarationID> >   SDEDeclaration;
        rule<ScannerT, parser_context<>, parser_tag<SDEModelDefID> >   SDEModelDef;
        rule<ScannerT, parser_context<>, parser_tag<modefullID> >         modefull;
        rule<ScannerT, parser_context<>, parser_tag<modesimplID> >        modesimpl;
        rule<ScannerT, parser_context<>, parser_tag<modeInterestID> >     modeInterest;

        rule<ScannerT, parser_context<>, parser_tag<ScriptID> > const&
        start() const { return Script; }
    };
  };
}

#endif

